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

#include "headers/pupildetection.h"

double calculateMean(const std::vector<double>& v)
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

std::vector<int> findGradientDirection(AOIProperties mAOI, double pupilXCentre, double pupilYCentre)
{
    std::vector<int> directionVector(mAOI.wdth * mAOI.hght);

    for (int y = 0; y < mAOI.hght; y++)
    {
        for (int x = 0; x < mAOI.wdth; x++)
        {
            double dx = x - pupilXCentre;
            double dy = pupilYCentre - y;

            int dir = calculateDirection(dx, dy);

            int centreIndex = mAOI.wdth * y + x;
            directionVector[centreIndex] = dir; // direction index of gradient
        }
    }

    return directionVector;
}

std::vector<int> radialGradient(const cv::Mat& img, const std::vector<int>& gradientDirections, int kernelSize)
{
    std::vector<int> dX = {  1,  1,  0, -1, -1, -1,  0,  1};
    std::vector<int> dY = {  0,  1,  1,  1,  0, -1, -1, -1};

    uchar *ptr = img.data;
    int width  = img.cols;
    int height = img.rows;

    std::vector<int> gradientVector(width * height, 0.0);

    int kernelRadius = (kernelSize - 1) / 2;

    for (int y = kernelRadius; y < height - kernelRadius; y++)
    {
        for (int x = kernelRadius; x < width - kernelRadius; x++)
        {
            int centreIndex = width * y + x;

            int direction = gradientDirections[centreIndex]; // direction index of gradient

            int dx = dX[direction];
            int dy = dY[direction];

            int neighbourIndexPositive = width * (y + dy * kernelRadius) + (x + dx * kernelRadius);
            int neighbourIndexNegative = width * (y - dy * kernelRadius) + (x - dx * kernelRadius);

            double positiveIntensity = ptr[neighbourIndexPositive];
            double negativeIntensity = ptr[neighbourIndexNegative];

            int val = positiveIntensity - negativeIntensity;
            if (val < 0) { val = 0; }

            gradientVector[centreIndex] = val;
        }
    }

    return gradientVector;
}

std::vector<int> nonMaximumSuppresion(const std::vector<int>& gradient, const std::vector<int>& gradientDirections, AOIProperties mAOI, double thresholdLow)
{
    std::vector<int> gradientSuppressed = gradient;

    std::vector<int> dX = {  1,  1,  0, -1, -1, -1,  0,  1};
    std::vector<int> dY = {  0,  1,  1,  1,  0, -1, -1, -1};

    for (int y = 0; y < mAOI.hght; y++)
    {
        for (int x = 0; x < mAOI.wdth; x++)
        {
            int centreIndex = mAOI.wdth * y + x;

            int gradientValue = gradient[centreIndex];

            if (gradientValue >= thresholdLow)
            {
                int direction = gradientDirections[centreIndex];

                int dx = dX[direction];
                int dy = dY[direction];

                int neighbourIndexPositive = mAOI.wdth * (y + dy) + (x + dx);
                int neighbourIndexNegative = mAOI.wdth * (y - dy) + (x - dx);

                if (gradientValue < gradient[neighbourIndexPositive] || gradientValue < gradient[neighbourIndexNegative])
                {
                    gradientSuppressed[centreIndex] = 0;
                }
            }
        }
    }

    return gradientSuppressed;
}

std::vector<int> hysteresisTracking(const std::vector<int>& gradientSuppressed, AOIProperties mAOI, int thresholdHigh, int thresholdLow)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

    int imgSize = mAOI.wdth * mAOI.hght;
    std::vector<int> edgePointVector(imgSize, 0);

    for (int iPixel = 0; iPixel < imgSize; iPixel++)
    {
        if (gradientSuppressed[iPixel] >= thresholdHigh && edgePointVector[iPixel] == 0)
        {
            edgePointVector[iPixel] = 1;
            std::vector<int> oldEdgeIndices = {iPixel};

            do
            {
                std::vector<int> newEdgeIndices;
                int numIndices = oldEdgeIndices.size();

                for (int iEdgePoint = 0; iEdgePoint < numIndices; iEdgePoint++)
                {
                    int centreIndex = oldEdgeIndices[iEdgePoint];
                    int centreXPos  = centreIndex % mAOI.wdth;
                    int centreYPos  = (centreIndex - centreXPos) / mAOI.wdth;

                    for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
                    {
                        int neighbourXPos = centreXPos + dX[m];
                        int neighbourYPos = centreYPos + dY[m];

                        if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght)
                        { continue; } // neighbour is out-of-bounds

                        int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

                        if (gradientSuppressed[neighbourIndex] >= thresholdLow && edgePointVector[neighbourIndex] != 1)
                        {
                            edgePointVector[neighbourIndex] = 1;
                            newEdgeIndices.push_back(neighbourIndex);
                        }
                    }
                }

                oldEdgeIndices = newEdgeIndices;

            } while (oldEdgeIndices.size() > 0);
        }
    }

    return edgePointVector;
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
                    int numberOfFilledPixels = 0;

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
                            numberOfFilledPixels++;
                        }
                    }

                    if (numberOfFilledPixels == 2) // if two neighbouring pixels in 4-connected environment are filled ...
                    {
                        int q = 2 * m + 1;

                        int xOpposite = xCentre + dX[q];
                        int yOpposite = yCentre + dY[q];
                        int iOpposite = mAOI.wdth * yOpposite + xOpposite;

                        if
                                (xOpposite < 0 || xOpposite >= mAOI.wdth  || // ... AND opposite pixel is out-of-bounds
                                 yOpposite < 0 || yOpposite >= mAOI.hght ||
                                 binaryImageVector[iOpposite] ==  0       || // ... OR unfilled ...
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

std::vector<int> getEdgeIndices(const std::vector<int>& binaryImageVector)
{
    int AOIArea = binaryImageVector.size();
    std::vector<int> cannyEdgeIndices;
    for (int iEdgePoint = 0; iEdgePoint < AOIArea; iEdgePoint++)
    { if (binaryImageVector[iEdgePoint] == 1) { cannyEdgeIndices.push_back(iEdgePoint); }}
    return cannyEdgeIndices;
}

inline bool findEdge(std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int& startEdgeIndex, double pupilXCentre, double pupilYCentre, bool USE_PRIOR_INFORMATION)
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

inline int connectEdges(std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int startIndex)
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
                edgePointNew = neighbourIndex; // edge point to-be-checked
                edgePoints.push_back(edgePointNew);
                break;
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
                int edgePointXPos = centreXPos + dX2[k];
                int edgePointYPos = centreYPos + dY2[k];
                if (edgePointXPos < 0 || edgePointXPos >= mAOI.wdth || edgePointYPos < 0 || edgePointYPos >= mAOI.hght) { continue; }
                int edgePointIndex = mAOI.wdth * edgePointYPos + edgePointXPos;

                int pointValue = cannyEdgeVector[edgePointIndex];

                if (pointValue == 1 || pointValue == 2)
                {
                    cannyEdgeVector[neighbourIndex] = 2;
                    return edgePointIndex;
                }
            }
        }
    }

    return startIndex;
}

inline void findEdgePoints(std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int startIndex)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

    cannyEdgeVector[startIndex] = 2; // tag pixel

    std::vector<int> edgePointsOld = {startIndex};

    int numEdgePoints;

    do
    {
        std::vector<int> edgePointsNew;

        numEdgePoints = edgePointsOld.size();

        for (int iEdgePoint = 0; iEdgePoint < numEdgePoints; iEdgePoint++) // loop through all newly added unchecked edge points
        {
            int centreIndex = edgePointsOld[iEdgePoint]; // index of current edge point

            int centreXPos = centreIndex % mAOI.wdth;
            int centreYPos = (centreIndex - centreXPos) / mAOI.wdth;

            bool NEW_EDGE_POINT_FOUND = false;

            for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
            {
                int neighbourXPos = centreXPos + dX[m];
                int neighbourYPos = centreYPos + dY[m];

                if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

                int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

                if (cannyEdgeVector[neighbourIndex] == 1) // if neighbouring point is filled ...
                {
                    cannyEdgeVector[neighbourIndex] = 2; // ... then tag it
                    edgePointsNew.push_back(neighbourIndex); // edge points to-be-checked
                    NEW_EDGE_POINT_FOUND = true;
                }
            }

            if (!NEW_EDGE_POINT_FOUND || iEdgePoint == 0)
            {
                int edgePointNew = connectEdges(cannyEdgeVector, mAOI, centreIndex); // Add additional point to the start or end of edge
                if (centreIndex != edgePointNew) { edgePointsNew.push_back(edgePointNew); }
            }
        }

        edgePointsOld = edgePointsNew;
        numEdgePoints  = edgePointsOld.size();
        edgePointsNew.clear();

    } while (numEdgePoints > 0);
}

inline void constructGraphTree(std::vector<vertexProperties>& verticesAll, std::vector<branchProperties>& branchesAll, std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int pointIndexStart)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

    vertexProperties vertexStart; // first vertex
    vertexStart.pointIndex = pointIndexStart;
    vertexStart.index      = 0;
    std::vector<vertexProperties> verticesOld = {vertexStart};

    int branchNumber = 0;
    int vertexNumber = 0;

    int numVertices;

    do
    {
        std::vector<vertexProperties> verticesNew; // store new found vertices

        numVertices = verticesOld.size();

        for (int iVertex = 0; iVertex < numVertices; iVertex++)
        {
            vertexProperties vertexCurrent = verticesOld[iVertex]; // run through all found vertices

            // Find all branches connected to current vertex

            int centreXPos =  vertexCurrent.pointIndex % mAOI.wdth;
            int centreYPos = (vertexCurrent.pointIndex - centreXPos) / mAOI.wdth;

            std::vector<int> edgePointIndicesStart;

            for (int m = 0; m < 8; m++) // loop through 8-connected environment of the vertex
            {
                int neighbourXPos = centreXPos + dX[m];
                int neighbourYPos = centreYPos + dY[m];

                if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

                int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

                if (cannyEdgeVector[neighbourIndex] == 2) // if neighbouring point is filled ...
                {
                    cannyEdgeVector[neighbourIndex] = 3; // ... then tag it
                    edgePointIndicesStart.push_back(neighbourIndex);
                }
            }

            // Loop through all found branches (if any) to find all edge points belonging to each branch

            int numEdges = edgePointIndicesStart.size();

            for (int iEdge = 0; iEdge < numEdges; iEdge++)
            {
                branchProperties branchNew;

                int edgePointIndexOld = edgePointIndicesStart[iEdge];

                while(true) // keep going until vertex is encountered
                {
                    std::vector<int> edgePointIndicesNew;

                    int centreIndex = edgePointIndexOld;
                    int centreXPos  = centreIndex % mAOI.wdth;
                    int centreYPos  = (centreIndex - centreXPos) / mAOI.wdth;

                    int neighbourIndex;

                    for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
                    {
                        int neighbourXPos = centreXPos + dX[m];
                        int neighbourYPos = centreYPos + dY[m];

                        if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds

                        neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;

                        if (cannyEdgeVector[neighbourIndex] == 2) // if neighbouring point was tagged previously ...
                        {
                            edgePointIndicesNew.push_back(neighbourIndex);
                        }
                    }

                    int numEdgePoints = edgePointIndicesNew.size();

                    if (numEdgePoints == 1) // edge continues
                    {
                        cannyEdgeVector[neighbourIndex] = 3; // ... give it a new tag
                        edgePointIndexOld = edgePointIndicesNew[0];
                        branchNew.pointIndices.push_back(edgePointIndexOld);
                        edgePointIndicesNew.clear();
                    }
                    else // new vertex found
                    {
                        vertexProperties vertexNew;
                        vertexNew.connections.push_back(branchNumber); // add current branch to new vertex connections
                        verticesNew.push_back(vertexNew); // new vertices to be checked
                        break;
                    }
                }

                vertexCurrent.connections.push_back(branchNumber); // add new found branches

                branchNew.length = branchNew.pointIndices.size();
                branchNew.index  = branchNumber;
                branchesAll.push_back(branchNew); // record all branches
                branchNumber++;
            }

            vertexCurrent.index = vertexNumber;
            verticesAll.push_back(vertexCurrent); // record all vertices
            vertexNumber++;
        }

        numVertices = verticesNew.size();
        verticesOld = verticesNew;

    } while (numVertices > 0);

    // Find vertices that branches are connected to

    numVertices = verticesAll.size();

    for (int iVertex = 0; iVertex < numVertices; iVertex++)
    {
        std::vector<int> vertexConnections = verticesAll[iVertex].connections;

        int numConnections = vertexConnections.size();

        for (int iEdge = 0; iEdge < numConnections; iEdge++)
        {
            branchesAll[vertexConnections[iEdge]].connections.push_back(iVertex);
        }
    }
}

inline std::vector<int> maxPath(const std::vector<std::vector<int>>& pathAll, const std::vector<branchProperties>& branchesAll)
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
            lengthTotal += branchesAll[pathCurrent[iBranch]].length;
        }

        pathLengths[iPath] = lengthTotal;
    }

    int pathIndexMax = std::distance(pathLengths.begin(), std::max_element(pathLengths.begin(), pathLengths.end()));
    return pathAll[pathIndexMax];
}

void findLongestPath(std::vector<int>& allPoints, std::vector<int>& pathPoints, std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int edgePointStart)
{
    std::vector<branchProperties> branchesAll;
    std::vector<vertexProperties> verticesAll;

    // Find all vertices and all connected branches (i.e. obtain graph tree)

    constructGraphTree(verticesAll, branchesAll, cannyEdgeVector, mAOI, edgePointStart);

    int numBranchesAll = branchesAll.size();
    int numVerticesAll = verticesAll.size();

    // Find longest path in edge collection

    std::vector<std::vector<int>> pathAll(numVerticesAll);

    for (int iVertex = 0; iVertex < numVerticesAll; iVertex++) // loop through all vertices
    {
        vertexProperties vertexStart = verticesAll[iVertex]; // first vertex

        std::vector<std::vector<int>> pathAllCurrent; // record all paths for current starting vertex

        std::vector<bool> branchFlagsAll(numBranchesAll,true);  // record if branches should be included again. Don't include terminals more than once

        while(true) // check all paths for current starting vertex
        {
            vertexProperties vertexNew = vertexStart;

            std::vector<bool> branchFlags(numBranchesAll,true); // no branch should be included twice

            std::vector<int> pathNew; // new path for current vertex

            bool PATH_CONTINUES = true;

            while(PATH_CONTINUES) // while path continues
            {
                int numConnections = vertexNew.connections.size();

                int branchNumber;
                PATH_CONTINUES = false;

                for (int iEdge = 0; iEdge < numConnections; iEdge++) // find first new branch that should still be included in path
                {
                    branchNumber = vertexNew.connections[iEdge];

                    if (branchFlagsAll[branchNumber] && branchFlags[branchNumber]) // branch should still be included
                    {
                        branchFlags[branchNumber] = false;
                        pathNew.push_back(branchNumber); // record branch indices for current path
                        PATH_CONTINUES = true;
                        break;
                    }
                }

                if (PATH_CONTINUES)
                {
                    // New vertex

                    std::vector<int> connectedVertices = branchesAll[branchNumber].connections; // get the two vertices the branch is connected to

                    for (int jVertex = 0; jVertex < 2; jVertex++) // find next vertex
                    {
                        int vertexIndex = connectedVertices[jVertex];

                        if (vertexIndex != vertexNew.index)
                        {
                            vertexNew = verticesAll[vertexIndex];
                            break;
                        }
                    }
                }
            }

            if (pathNew.size() > 0)
            {
                branchFlagsAll[pathNew.back()] = false; // don't include last included branch again
                pathAllCurrent.push_back(pathNew);
            }
            else { break; } // no new path found
        }

        // Find longest path for current starting vertex

        pathAll[iVertex] = maxPath(pathAllCurrent,branchesAll); // record longest path
    }

    // Longest path overall

    std::vector<int> pathBranchIndices = maxPath(pathAll,branchesAll);
    int numBranchesPath = pathBranchIndices.size();

    // Combine all edge point indices of longest path

    for (int iBranch = 0; iBranch < numBranchesPath; iBranch++)
    {
        int branchIndex = pathBranchIndices[iBranch];
        pathPoints.insert(std::end(pathPoints), std::begin(branchesAll[branchIndex].pointIndices), std::end(branchesAll[branchIndex].pointIndices));
    }

    // All edge points

    for (int iBranch = 0; iBranch < numBranchesAll; iBranch++)
    {
        allPoints.insert(std::end(allPoints), std::begin(branchesAll[iBranch].pointIndices), std::end(branchesAll[iBranch].pointIndices));
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
        findEdgePoints(cannyEdgeVector, mAOI, startIndex); // tag all desired edges
        startIndices.push_back(startIndex);
    }

    int numEdges = startIndices.size();

    for (int iEdge = 0; iEdge < numEdges; iEdge++)
    {
        // Find longest path

        std::vector<int> pathIndices;
        std::vector<int>  allIndices;

        findLongestPath(allIndices, pathIndices, cannyEdgeVector, mAOI, startIndices[iEdge]);

        // Give pixels that have been added to new outline a new tag

        for (int iEdgePoint = 0, edgeLength = pathIndices.size(); iEdgePoint < edgeLength; iEdgePoint++)
        { cannyEdgeVector[pathIndices[iEdgePoint]] = 4; }

        // Remove tag from pixels that have been tagged before, but not included in new outline

        for (int iEdgePoint = 0, edgeLength = allIndices.size(); iEdgePoint < edgeLength; iEdgePoint++)
        {
            int edgePointIndex = allIndices[iEdgePoint];

            if (cannyEdgeVector[edgePointIndex] == 2 || cannyEdgeVector[edgePointIndex] == 3)
            {  cannyEdgeVector[edgePointIndex] = 1; }
        }

        int edgeLength = pathIndices.size();

        edgeProperties mEdgeProperties;
        mEdgeProperties.length       = edgeLength;
        mEdgeProperties.pointIndices = pathIndices;

        vEdgePropertiesAll.push_back(mEdgeProperties);
    }

    return vEdgePropertiesAll;
}

std::vector<double> calculateCurvatures(std::vector<double>& edgeXNormals, std::vector<double>& edgeYNormals, const std::vector<double>& edgeXTangents, const std::vector<double>& edgeYTangents)
{
    int numPos = 0;
    int numNeg = 0;

    std::vector<double> edgeCurvatures;

    int edgeLength = edgeXNormals.size();

    for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeLength - curvatureWindowLength; iEdgePoint++)
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

        edgeCurvatures[iEdgePoint] = 180 * vectorAngle / M_PI; // in degrees

        edgeXNormals[iEdgePoint] = -meanXTangent_1 + meanXTangent_2;
        edgeYNormals[iEdgePoint] = -meanYTangent_1 + meanYTangent_2;
    }

    if (numNeg > numPos) // if majority sign is negative, then swap all signs
    {
        for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeLength - curvatureWindowLength; iEdgePoint++)
        {
            edgeCurvatures[iEdgePoint] = -edgeCurvatures[iEdgePoint];
        }
    }

    return edgeCurvatures;
}

std::vector<edgeProperties> edgeSegmentation(edgeProperties mEdgeProperties, double curvatureLowerLimit, double curvatureUpperLimit)
{
    std::vector<int> edgePointsAll      = mEdgeProperties.pointIndices;
    std::vector<double> edgeCurvatures  = mEdgeProperties.curvatures;
    std::vector<double> edgeIntensities = mEdgeProperties.intensities;

    int edgeLength = edgePointsAll.size();

    std::vector<edgeProperties> vEdgePropertiesNew;

    // find breakpoints

    std::vector<int> breakPointIndices; // position of breakpoints
    breakPointIndices.push_back(0); // add first point

    for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeLength - curvatureWindowLength; iEdgePoint++)
    {
        double curvature = edgeCurvatures[iEdgePoint];

        if (std::abs(curvature) >= curvatureUpperLimit || curvature <= curvatureLowerLimit)
        {
            breakPointIndices.push_back(iEdgePoint);
        }
    }

    breakPointIndices.push_back(edgeLength - 1); // add last point

    // evaluate each partial edge

    std::sort(breakPointIndices.begin(), breakPointIndices.end());

    int numBreakPoints = breakPointIndices.size();

    for (int iBreakPoint = 0; iBreakPoint < numBreakPoints - 1; iBreakPoint++)
    {
        int iStartBreakPoint = breakPointIndices[iBreakPoint] + 1;
        int iEndBreakPoint   = breakPointIndices[iBreakPoint  + 1];
        int edgeLengthNew    = iEndBreakPoint - iStartBreakPoint;

        // grab indices of edge points

        if (edgeLengthNew < curvatureWindowLength) { continue; } // ignore short edges
        std::vector<int>    edgePointsNew     (edgePointsAll.begin()   + iStartBreakPoint, edgePointsAll.begin()   + iEndBreakPoint);
        std::vector<double> edgeCurvaturesNew (edgeCurvatures.begin()  + iStartBreakPoint, edgeCurvatures.begin()  + iEndBreakPoint);
        std::vector<double> edgeIntensitiesNew(edgeIntensities.begin() + iStartBreakPoint, edgeIntensities.begin() + iEndBreakPoint);

        edgeProperties mEdgeProperties;
        mEdgeProperties.length       = edgeLengthNew;
        mEdgeProperties.pointIndices = edgePointsNew;
        mEdgeProperties.curvatures   = edgeCurvaturesNew;
        mEdgeProperties.intensities  = edgeIntensitiesNew;

        vEdgePropertiesNew.push_back(mEdgeProperties);
    }

    return vEdgePropertiesNew;
}

std::vector<double> findEdgeIntensities(const cv::Mat& img, const std::vector<int>& edgePoints, const std::vector<double>& edgePointXNormals, const std::vector<double>& edgePointYNormals,  AOIProperties mAOI)
{
    int edgeLength = edgePoints.size();

    // calculate pixel intensities within inner curve of edge points

    std::vector<double> edgeIntensities(edgeLength);

    uchar *ptr_img = img.data;

    for (int iEdgePoint = 0; iEdgePoint < edgeLength; iEdgePoint++)
    {
        int edgePointIndex = edgePoints[iEdgePoint];
        int edgePointXPos  = edgePointIndex % mAOI.wdth;
        int edgePointYPos  = (edgePointIndex - edgePointXPos) / mAOI.wdth;

        int offsetXPos = edgePointXPos + edgeIntensitiesPositionOffset * ceil2(edgePointXNormals[iEdgePoint]);
        int offsetYPos = edgePointYPos + edgeIntensitiesPositionOffset * ceil2(edgePointYNormals[iEdgePoint]);

        if (offsetXPos < 0 || offsetXPos > mAOI.wdth || offsetYPos < 0 || offsetYPos > mAOI.hght)
        {       edgeIntensities[iEdgePoint] = ptr_img[edgePointXPos + edgePointYPos * mAOI.wdth]; }
        else {  edgeIntensities[iEdgePoint] = ptr_img[   offsetXPos +    offsetYPos * mAOI.wdth]; }
    }

    return edgeIntensities;
}

void calculateCurvatureRange(edgeProperties& mEdgeProperties)
{
    double curvatureAvg = 0;
    double curvatureMax;
    double curvatureMin;

    std::vector<double> edgeCurvaturesNew;

    for (int iEdgePoint = 0; iEdgePoint < mEdgeProperties.length; iEdgePoint++)
    {
        double curvature = mEdgeProperties.curvatures[iEdgePoint];
        if (curvature < 180.0) { edgeCurvaturesNew.push_back(curvature); }
    }

    int edgeLengthNew = edgeCurvaturesNew.size();

    if (edgeLengthNew > 0)
    {
        for (int iEdgePoint = 0; iEdgePoint < edgeLengthNew; iEdgePoint++)
        {
            curvatureAvg += std::abs(edgeCurvaturesNew[iEdgePoint] / edgeLengthNew);
        }

        curvatureMax = *std::max_element(std::begin(edgeCurvaturesNew), std::end(edgeCurvaturesNew));
        curvatureMin = *std::min_element(std::begin(edgeCurvaturesNew), std::end(edgeCurvaturesNew));
    }
    else
    {
        curvatureAvg = 360;
        curvatureMax = 360;
        curvatureMin = 360;
    }

    mEdgeProperties.curvatureAvg = curvatureAvg;
    mEdgeProperties.curvatureMax = curvatureMax;
    mEdgeProperties.curvatureMin = curvatureMin;
}

void restoreEdgePoints(edgeProperties& mEdgeProperties, std::vector<int>& cannyEdgeVector, AOIProperties mAOI)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};

    int edgeLength = mEdgeProperties.length;

    // add additional adjacent indices that were removed by morphological operation

    for (int iEdgePoint = 0; iEdgePoint < edgeLength; iEdgePoint++)
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
                cannyEdgeVector[neighbourIndex] = 4; // ... tag it and ...
                mEdgeProperties.pointIndices.push_back(neighbourIndex); // ... add it to the (partial) edge
            }
        }
    }

    mEdgeProperties.size = mEdgeProperties.pointIndices.size();
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
        double dCurvature       = mDetectionProperties.v.edgeCurvaturePrediction - vEdgePropertiesAll[iEdge].curvatureAvg;
        double intensity        = vEdgePropertiesAll[iEdge].intensityAvg;

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
        std::sort(totalScoresSorted.begin(), totalScoresSorted.end());
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
    mEllipseProperties.pupilDetected = true;
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
        if (std::isnan(ellipseParameters[iParameter])) { mEllipseProperties.pupilDetected = false; } // ERROR
    }

    return mEllipseProperties;
}

ellipseProperties findBestEllipseFit(const std::vector<edgeProperties>& vEdgePropertiesAll, AOIProperties mAOI, detectionProperties mDetectionProperties, bool USE_PRIOR_INFORMATION)
{
    ellipseProperties mEllipseProperties;
    mEllipseProperties.pupilDetected = false;

    int totalNumberOfEdges = vEdgePropertiesAll.size(); // total number of edges

    std::vector<ellipseProperties> vEllipsePropertiesAll; // vector to record information for each accepted ellipse fit

    int numFits = 0; // iterator over all accepted combinations

    for (int combiNumEdges = totalNumberOfEdges; combiNumEdges >= 1; combiNumEdges--) // loop through all possible edge set sizes
    {
        std::vector<bool> edgeCombination(totalNumberOfEdges);
        std::fill(edgeCombination.begin() + totalNumberOfEdges - combiNumEdges, edgeCombination.end(), true);

        do // loop through all possible edge combinations for the current set size
        {
            std::vector<double> combiEdgeIntensities(combiNumEdges);
            std::vector<int>    combiEdgeIndices    (combiNumEdges);
            std::vector<int>    combiEdgeLengths    (combiNumEdges);
            std::vector<int>    combiEdgeSizes      (combiNumEdges);
            std::vector<std::vector<int>> combiEdgePointIndices(combiNumEdges);

            for (int iEdge = 0, jEdge = 0; iEdge < totalNumberOfEdges; ++iEdge)
            {
                if (edgeCombination[iEdge])
                {
                    combiEdgeIndices     [jEdge] = vEdgePropertiesAll[iEdge].index;
                    combiEdgeIntensities [jEdge] = vEdgePropertiesAll[iEdge].intensityAvg;
                    combiEdgeLengths     [jEdge] = vEdgePropertiesAll[iEdge].length;
                    combiEdgeSizes       [jEdge] = vEdgePropertiesAll[iEdge].size;
                    combiEdgePointIndices[jEdge] = vEdgePropertiesAll[iEdge].pointIndices;
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
            { edgeIndices.insert(edgeIndices.end(), combiEdgePointIndices[iEdge].begin(), combiEdgePointIndices[iEdge].end()); }

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

            if (!mEllipsePropertiesNew.pupilDetected) { continue; } // error

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
            mEllipsePropertiesNew.intensity   = calculateMean(combiEdgeIntensities);
            mEllipsePropertiesNew.edgeLength  = edgeSetLength;
            vEllipsePropertiesAll.push_back(mEllipsePropertiesNew);

            numFits++;
        }
        while (std::next_permutation(edgeCombination.begin(), edgeCombination.end()));
    }

    // grab ellipse fit that resembles prior the most in size and shape

    if (numFits > 0)
    {
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
            double dx = vEllipsePropertiesAll[iFit].xPos - mDetectionProperties.v.xPosPrediction;
            double dy = vEllipsePropertiesAll[iFit].yPos - mDetectionProperties.v.yPosPrediction;
            double displacementChange  = sqrt(pow(dx,2) + pow(dy,2));
            double circumferenceChange = (std::abs(vEllipsePropertiesAll[iFit].circumference - mDetectionProperties.v.circumferencePrediction));
            double aspectRatioChange   = (std::abs(vEllipsePropertiesAll[iFit].aspectRatio   - mDetectionProperties.v.aspectRatioPrediction));
            double lengthDifference    = (std::abs(vEllipsePropertiesAll[iFit].edgeLength    - mDetectionProperties.v.circumferencePrediction));
            double fitError            = vEllipsePropertiesAll[iFit].fitError;

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

        mEllipseProperties = vEllipsePropertiesAll[acceptedFitIndex];
        mEllipseProperties.pupilDetected = true;
    }

    return mEllipseProperties;
}

detectionProperties pupilDetection(const cv::Mat& imageOriginalBGR, detectionProperties mDetectionProperties)
{
    // Define some variables

    detectionProperties mDetectionPropertiesNew = mDetectionProperties; // new properties for new frame
    std::vector<edgeProperties> vEdgePropertiesNew;
    ellipseProperties mEllipseProperties;

    mEllipseProperties.pupilDetected        = false;
    mDetectionPropertiesNew.m.errorDetected = false;

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

    AOIProperties outerAOI;
    outerAOI.xPos = 0;
    outerAOI.yPos = 0;
    outerAOI.wdth = 0;
    outerAOI.hght = 0;

    if (innerAOI.wdth > 0 && innerAOI.hght > 0)
    {
        // Needed for offline mode when threshold is too low

        if (mDetectionProperties.v.thresholdDisplacementChange > imageWdth)
        {   mDetectionProperties.v.thresholdDisplacementChange = imageWdth; }
        else if (mDetectionProperties.v.thresholdDisplacementChange < mDetectionProperties.p.displacementChangeThreshold)
        {        mDetectionProperties.v.thresholdDisplacementChange = mDetectionProperties.p.displacementChangeThreshold; }

        if (mDetectionProperties.v.thresholdCircumferenceChange > mDetectionProperties.p.circumferenceMax)
        {   mDetectionProperties.v.thresholdCircumferenceChange = mDetectionProperties.p.circumferenceMax; }
        else if (mDetectionProperties.v.thresholdCircumferenceChange < mDetectionProperties.p.circumferenceChangeThreshold)
        {        mDetectionProperties.v.thresholdCircumferenceChange = mDetectionProperties.p.circumferenceChangeThreshold; }

        if (mDetectionProperties.v.thresholdAspectRatioChange > 1.0)
        {   mDetectionProperties.v.thresholdAspectRatioChange = 1.0; }
        else if (mDetectionProperties.v.thresholdAspectRatioChange < mDetectionProperties.p.aspectRatioChangeThreshold)
        {        mDetectionProperties.v.thresholdAspectRatioChange = mDetectionProperties.p.aspectRatioChangeThreshold; }

        if (mDetectionProperties.v.curvatureOffset > 180)
        {   mDetectionProperties.v.curvatureOffset = 180; }
        else if (mDetectionProperties.v.curvatureOffset < mDetectionProperties.p.curvatureOffset)
        {        mDetectionProperties.v.curvatureOffset = mDetectionProperties.p.curvatureOffset; }

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

        double xPosPredictionRelative = mDetectionProperties.v.xPosPrediction - outerAOI.xPos;
        double yPosPredictionRelative = mDetectionProperties.v.yPosPrediction - outerAOI.yPos;

        std::vector<int> cannyEdges; // binary vector

        if (USE_PRIOR_INFORMATION)
        {
            std::vector<int>  dirGradient           = findGradientDirection(outerAOI, xPosPredictionRelative, yPosPredictionRelative);
            std::vector<int>  imgGradient           = radialGradient(imageAOIGrayBlurred, dirGradient, mDetectionProperties.p.cannyKernelSize);
            std::vector<int>  imgGradientSuppressed = nonMaximumSuppresion(imgGradient, dirGradient, outerAOI, mDetectionProperties.p.cannyThresholdLow);
            cannyEdges = hysteresisTracking(imgGradientSuppressed, outerAOI, mDetectionProperties.p.cannyThresholdHigh, mDetectionProperties.p.cannyThresholdLow);
        }
        else
        {
            cv::Mat imageCannyEdges;
            cv::Canny(imageAOIGrayBlurred, imageCannyEdges, 30 * mDetectionProperties.p.cannyThresholdHigh, 30 * mDetectionProperties.p.cannyThresholdLow, 5);
            cannyEdges = cannyConversion(imageCannyEdges, outerAOI);
        }

        std::vector<int> cannyEdgesSharpened = sharpenEdges(cannyEdges, outerAOI); // Morphological operation
        std::vector<int> edgeIndices         = getEdgeIndices(cannyEdges);

        //////////////////////////////////////////////////////////////////////////////
        /////////////////////// EDGE SELECTION & SEGMENTATION  ///////////////////////
        //////////////////////////////////////////////////////////////////////////////

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

        std::vector<edgeProperties> vEdgePropertiesAll = edgeSelection(cannyEdgesSharpened, outerAOI, xPosPredictionRelative, yPosPredictionRelative, USE_PRIOR_INFORMATION);

        int numEdges = vEdgePropertiesAll.size();
        std::vector<edgeProperties> vEdgePropertiesNew;

        for (int iEdge = 0; iEdge < numEdges; iEdge++) // ignore short edges
        {
            edgeProperties mEdgeProperties = vEdgePropertiesAll[iEdge];
            int edgeLength = mEdgeProperties.length;
            if (edgeLength > curvatureWindowLength) { vEdgePropertiesNew.push_back(mEdgeProperties); }
        }

        vEdgePropertiesAll.clear();

        int numEdgesNew = vEdgePropertiesNew.size();

        for (int iEdge = 0; iEdge < numEdgesNew; iEdge++)
        {
            edgeProperties mEdgeProperties = vEdgePropertiesNew[iEdge];

            int edgeLength = mEdgeProperties.length;

            std::vector<double> edgeXTangents(edgeLength);
            std::vector<double> edgeYTangents(edgeLength);

            calculateEdgeDirections(mEdgeProperties.pointIndices, edgeXTangents, edgeYTangents, outerAOI);

            std::vector<double> edgePointXNormals(edgeLength);
            std::vector<double> edgePointYNormals(edgeLength);

            mEdgeProperties.curvatures  = calculateCurvatures(edgePointXNormals, edgePointYNormals, edgeXTangents, edgeYTangents);
            mEdgeProperties.intensities = findEdgeIntensities(imageAOIGray, mEdgeProperties.pointIndices, edgePointXNormals, edgePointYNormals, outerAOI);

            std::vector<edgeProperties> vEdgePropertiesTemp = edgeSegmentation(mEdgeProperties, curvatureLowerLimit, curvatureUpperLimit);

            int numEdgesTemp = vEdgePropertiesTemp.size();
            for (int jEdge = 0; jEdge < numEdgesTemp; jEdge++) { restoreEdgePoints(vEdgePropertiesTemp[jEdge], cannyEdgesSharpened, outerAOI); }

            vEdgePropertiesAll.insert(vEdgePropertiesNew.end(), vEdgePropertiesTemp.begin(), vEdgePropertiesTemp.end()); // record edges
        }

        // Remove short edges

        ///////////////////////////////////////////////////////////////////////
        /////////////////////// EDGE CLASSIFICATION  //////////////////////////
        ///////////////////////////////////////////////////////////////////////

        int numEdgesTotal = vEdgePropertiesAll.size();

        if (numEdgesTotal > 0) // THRESHOLD: ignore empty edge collections
        {
            for (int iEdge = 0; iEdge < numEdgesTotal; iEdge++) // initialize edges
            {
                vEdgePropertiesAll[iEdge].flag  = 0;
                vEdgePropertiesAll[iEdge].index = iEdge;
            }

            for (int iEdge = 0; iEdge < numEdgesTotal; iEdge++) // calculate distance between each edge point and expected pupil centre
            {
                double dR = 0;
                int edgeSize = vEdgePropertiesAll[iEdge].size;

                for (int iEdgePoint = 0; iEdgePoint < edgeSize; iEdgePoint++)
                {
                    int edgePointIndex = vEdgePropertiesAll[iEdge].pointIndices[iEdgePoint];
                    int edgePointXPos  =  edgePointIndex % outerAOI.wdth;
                    int edgePointYPos  = (edgePointIndex - edgePointXPos) / outerAOI.wdth;

                    double dX = xPosPredictionRelative - edgePointXPos;
                    double dY = yPosPredictionRelative - edgePointYPos;

                    dR += sqrt(pow(dX, 2) + pow(dY, 2));
                }

                vEdgePropertiesAll[iEdge].radius = dR / edgeSize; // average radius
            }

            mDetectionProperties.v.edgeCurvaturePrediction = 0.5 * (curvatureUpperLimit + curvatureLowerLimit);

            std::vector<int> acceptedEdges = edgeClassification(mDetectionProperties, vEdgePropertiesAll, USE_PRIOR_INFORMATION);
            int numEdges = acceptedEdges.size();

            for (int iEdge = 0; iEdge < numEdges; iEdge++) // grab accepted edges
            {
                int jEdge = acceptedEdges[iEdge];
                vEdgePropertiesAll[jEdge].flag = 1; // tag it
                vEdgePropertiesNew.push_back(vEdgePropertiesAll[jEdge]);
            }

            //////////////////////////////////////////////////////////////////
            /////////////////////// ELLIPSE FITTING  /////////////////////////
            //////////////////////////////////////////////////////////////////

            auto t1 = std::chrono::high_resolution_clock::now();
            mEllipseProperties = findBestEllipseFit(vEdgePropertiesNew, outerAOI, mDetectionProperties, USE_PRIOR_INFORMATION); // ellipse fitting
            auto t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
            std::cout << fp_ms.count() << std::endl;

            // Classify edges

            int numEdgesNew = mEllipseProperties.edgeIndices.size();

            for (int iEdge = 0; iEdge < numEdgesNew; iEdge++)
            {
                int jEdge = mEllipseProperties.edgeIndices[iEdge];
                vEdgePropertiesAll[jEdge].flag = 2;
            }
        }

        /////////////////////////////////////////////////////////////////
        /////////////////////// SAVING DATA  ////////////////////////////
        /////////////////////////////////////////////////////////////////

        // Features for all edges

        mDetectionPropertiesNew.m.edgePropertiesAll = vEdgePropertiesAll;

        // Save parameters

        mDetectionPropertiesNew.v.edgeCurvaturePrediction = mDetectionProperties.v.edgeCurvaturePrediction;

        mDetectionPropertiesNew.v.pupilDetected = mEllipseProperties.pupilDetected;

        // For draw functions

        mDetectionPropertiesNew.m.outerAOI = outerAOI;
        mDetectionPropertiesNew.m.innerAOI = innerAOI;
        mDetectionPropertiesNew.m.glintAOI = glintAOI;

        mDetectionPropertiesNew.m.cannyEdgeIndices    = edgeIndices;
        mDetectionPropertiesNew.m.ellipseCoefficients = mEllipseProperties.coefficients;
    }
    else
    {
        mDetectionPropertiesNew.m.errorDetected = true;
    }

    // For running averages

    if (!mEllipseProperties.pupilDetected) // pupil not detected
    {
        mDetectionPropertiesNew.v.aspectRatioAverage    = mDetectionProperties.v.aspectRatioAverage + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.aspectRatioPrediction - mDetectionProperties.v.aspectRatioAverage);
        mDetectionPropertiesNew.v.aspectRatioMomentum   = mDetectionProperties.v.aspectRatioMomentum * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.aspectRatioPrediction = mDetectionProperties.v.aspectRatioPrediction + mDetectionProperties.p.alphaPrediction * (mDetectionPropertiesNew.v.aspectRatioAverage - mDetectionProperties.v.aspectRatioPrediction);

        mDetectionPropertiesNew.v.circumferenceAverage    = mDetectionProperties.v.circumferenceAverage + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.circumferencePrediction - mDetectionProperties.v.circumferenceAverage);
        mDetectionPropertiesNew.v.circumferenceMomentum   = mDetectionProperties.v.circumferenceMomentum * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.circumferencePrediction = mDetectionProperties.v.circumferencePrediction + mDetectionProperties.p.alphaPrediction * (mDetectionPropertiesNew.v.circumferenceAverage - mDetectionProperties.v.circumferencePrediction);

        mDetectionPropertiesNew.v.curvatureOffset = mDetectionProperties.v.curvatureOffset * (1 / mDetectionProperties.p.alphaMiscellaneous);

        mDetectionPropertiesNew.v.heightAverage    = mDetectionProperties.v.heightAverage + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.heightPrediction - mDetectionProperties.v.heightAverage);
        mDetectionPropertiesNew.v.heightMomentum   = mDetectionProperties.v.heightMomentum * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.heightPrediction = mDetectionProperties.v.heightPrediction + mDetectionProperties.p.alphaPrediction * (mDetectionPropertiesNew.v.heightAverage - mDetectionProperties.v.heightPrediction);

        mDetectionPropertiesNew.v.radiusMomentum   = mDetectionProperties.v.radiusMomentum * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.radiusPrediction = mDetectionProperties.v.circumferencePrediction / (2 * M_PI);

        mDetectionPropertiesNew.v.searchRadius = mDetectionProperties.v.searchRadius * (1 / mDetectionProperties.p.alphaMiscellaneous);

        mDetectionPropertiesNew.v.thresholdAspectRatioChange   = mDetectionProperties.v.thresholdAspectRatioChange   * (1 / mDetectionProperties.p.alphaMiscellaneous);
        mDetectionPropertiesNew.v.thresholdCircumferenceChange = mDetectionProperties.v.thresholdCircumferenceChange * (1 / mDetectionProperties.p.alphaMiscellaneous);
        mDetectionPropertiesNew.v.thresholdDisplacementChange  = mDetectionProperties.v.thresholdDisplacementChange * (1 / mDetectionProperties.p.alphaMiscellaneous);

        mDetectionPropertiesNew.v.widthAverage    = mDetectionProperties.v.widthAverage + mDetectionProperties.p.alphaAverage * (mDetectionProperties.v.widthPrediction - mDetectionProperties.v.widthAverage);
        mDetectionPropertiesNew.v.widthMomentum   = mDetectionProperties.v.widthMomentum * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.widthPrediction = mDetectionProperties.v.widthPrediction + mDetectionProperties.p.alphaPrediction * (mDetectionPropertiesNew.v.widthAverage - mDetectionProperties.v.widthPrediction);

        mDetectionPropertiesNew.v.xPosPrediction = mDetectionProperties.v.xPosPrediction + mDetectionProperties.p.alphaPrediction * (outerAOI.xPos + 0.5 * outerAOI.wdth - mDetectionProperties.v.xPosPrediction) + mDetectionProperties.v.xVelocity;
        mDetectionPropertiesNew.v.yPosPrediction = mDetectionProperties.v.yPosPrediction + mDetectionProperties.p.alphaPrediction * (outerAOI.yPos + 0.5 * outerAOI.hght - mDetectionProperties.v.yPosPrediction) + mDetectionProperties.v.yVelocity;

        mDetectionPropertiesNew.v.xVelocity     = mDetectionProperties.v.xVelocity * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.yVelocity     = mDetectionProperties.v.yVelocity * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.priorCertainty = mDetectionProperties.v.priorCertainty * mDetectionProperties.p.alphaMiscellaneous;
    }
    else // pupil detected
    {
        mDetectionPropertiesNew.v.aspectRatioExact   = mEllipseProperties.aspectRatio;
        mDetectionPropertiesNew.v.circumferenceExact = mEllipseProperties.circumference;

        mDetectionPropertiesNew.v.aspectRatioPrediction =  mDetectionProperties.v.aspectRatioPrediction + mDetectionProperties.p.alphaPrediction * (mEllipseProperties.aspectRatio - mDetectionProperties.v.aspectRatioPrediction) + mDetectionProperties.v.aspectRatioMomentum;
        mDetectionPropertiesNew.v.aspectRatioAverage    =  mDetectionProperties.v.aspectRatioAverage    + mDetectionProperties.p.alphaAverage    * (mDetectionProperties.v.aspectRatioPrediction - mDetectionProperties.v.aspectRatioAverage);
        mDetectionPropertiesNew.v.aspectRatioMomentum   = (mDetectionProperties.v.aspectRatioMomentum   + mDetectionPropertiesNew.v.aspectRatioPrediction - mDetectionProperties.v.aspectRatioPrediction) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.circumferencePrediction =  mDetectionProperties.v.circumferencePrediction + mDetectionProperties.p.alphaPrediction * (mEllipseProperties.circumference - mDetectionProperties.v.circumferencePrediction) + mDetectionProperties.v.circumferenceMomentum;
        mDetectionPropertiesNew.v.circumferenceAverage    =  mDetectionProperties.v.circumferenceAverage    + mDetectionProperties.p.alphaAverage    * (mDetectionProperties.v.circumferencePrediction - mDetectionProperties.v.circumferenceAverage);
        mDetectionPropertiesNew.v.circumferenceMomentum   = (mDetectionProperties.v.circumferenceMomentum   + mDetectionPropertiesNew.v.circumferencePrediction - mDetectionProperties.v.circumferencePrediction) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.curvatureOffset = mDetectionProperties.v.curvatureOffset * mDetectionProperties.p.alphaMiscellaneous;

        mDetectionPropertiesNew.v.heightPrediction =  mDetectionProperties.v.heightPrediction + mDetectionProperties.p.alphaPrediction * (mEllipseProperties.height - mDetectionProperties.v.heightPrediction) + mDetectionProperties.v.heightMomentum;
        mDetectionPropertiesNew.v.heightAverage    =  mDetectionProperties.v.heightAverage    + mDetectionProperties.p.alphaAverage    * (mDetectionProperties.v.heightPrediction - mDetectionProperties.v.heightAverage);
        mDetectionPropertiesNew.v.heightMomentum   = (mDetectionProperties.v.heightMomentum   + mDetectionPropertiesNew.v.heightPrediction - mDetectionProperties.v.heightPrediction) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.radiusPrediction =  mDetectionProperties.v.radiusPrediction + mDetectionProperties.p.alphaPrediction * (mEllipseProperties.radius - mDetectionProperties.v.radiusPrediction) + mDetectionProperties.v.radiusMomentum;
        mDetectionPropertiesNew.v.radiusMomentum   = (mDetectionProperties.v.radiusMomentum + (mDetectionPropertiesNew.v.radiusPrediction - mDetectionProperties.v.radiusPrediction)) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.searchRadius = mDetectionProperties.v.searchRadius * mDetectionProperties.p.alphaMiscellaneous;

        mDetectionPropertiesNew.v.thresholdAspectRatioChange   = mDetectionProperties.v.thresholdAspectRatioChange   * mDetectionProperties.p.alphaMiscellaneous;
        mDetectionPropertiesNew.v.thresholdCircumferenceChange = mDetectionProperties.v.thresholdCircumferenceChange * mDetectionProperties.p.alphaMiscellaneous;
        mDetectionPropertiesNew.v.thresholdDisplacementChange  = mDetectionProperties.v.thresholdDisplacementChange  * mDetectionProperties.p.alphaMiscellaneous;

        mDetectionPropertiesNew.v.widthPrediction =  mDetectionProperties.v.widthPrediction + mDetectionProperties.p.alphaPrediction * (mEllipseProperties.width - mDetectionProperties.v.widthPrediction) + mDetectionProperties.v.widthMomentum;
        mDetectionPropertiesNew.v.widthAverage    =  mDetectionProperties.v.widthAverage    + mDetectionProperties.p.alphaAverage    * (mDetectionProperties.v.widthPrediction - mDetectionProperties.v.widthAverage);
        mDetectionPropertiesNew.v.widthMomentum   = (mDetectionProperties.v.widthMomentum   + mDetectionPropertiesNew.v.widthPrediction - mDetectionProperties.v.widthPrediction) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.xPosExact = mEllipseProperties.xPos + outerAOI.xPos;
        mDetectionPropertiesNew.v.yPosExact = mEllipseProperties.yPos + outerAOI.yPos;

        mDetectionPropertiesNew.v.xPosPrediction = mDetectionProperties.v.xPosPrediction + mDetectionProperties.p.alphaPrediction * (mDetectionPropertiesNew.v.xPosExact - mDetectionProperties.v.xPosPrediction) + mDetectionProperties.v.xVelocity;
        mDetectionPropertiesNew.v.yPosPrediction = mDetectionProperties.v.yPosPrediction + mDetectionProperties.p.alphaPrediction * (mDetectionPropertiesNew.v.yPosExact - mDetectionProperties.v.yPosPrediction) + mDetectionProperties.v.yVelocity;

        mDetectionPropertiesNew.v.xVelocity = (mDetectionProperties.v.xVelocity + mDetectionPropertiesNew.v.xPosPrediction - mDetectionProperties.v.xPosPrediction) * mDetectionProperties.p.alphaMomentum;
        mDetectionPropertiesNew.v.yVelocity = (mDetectionProperties.v.yVelocity + mDetectionPropertiesNew.v.yPosPrediction - mDetectionProperties.v.yPosPrediction) * mDetectionProperties.p.alphaMomentum;

        mDetectionPropertiesNew.v.priorCertainty = mDetectionProperties.v.priorCertainty * (1 / mDetectionProperties.p.alphaMiscellaneous);

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

    // Check variable limits

    if (mDetectionPropertiesNew.v.searchRadius > 0.5 * imageWdth)
    {   mDetectionPropertiesNew.v.searchRadius = 0.5 * imageWdth; }
    else if (mDetectionPropertiesNew.v.searchRadius < (0.5 * outerAOI.wdth))
    {        mDetectionPropertiesNew.v.searchRadius = ceil(0.5 * outerAOI.wdth); }

    if (mDetectionPropertiesNew.v.searchRadius > 0.5 * imageHght)
    {   mDetectionPropertiesNew.v.searchRadius = 0.5 * imageHght; }
    else if (mDetectionPropertiesNew.v.searchRadius < (0.5 * outerAOI.hght))
    {        mDetectionPropertiesNew.v.searchRadius = ceil(0.5 * outerAOI.hght); }

    if (mDetectionPropertiesNew.v.thresholdAspectRatioChange > 1.0)
    {   mDetectionPropertiesNew.v.thresholdAspectRatioChange = 1.0; }
    else if (mDetectionPropertiesNew.v.thresholdAspectRatioChange < mDetectionPropertiesNew.p.aspectRatioChangeThreshold)
    {        mDetectionPropertiesNew.v.thresholdAspectRatioChange = mDetectionPropertiesNew.p.aspectRatioChangeThreshold; }

    if (mDetectionPropertiesNew.v.thresholdCircumferenceChange > mDetectionPropertiesNew.p.circumferenceMax)
    {   mDetectionPropertiesNew.v.thresholdCircumferenceChange = mDetectionPropertiesNew.p.circumferenceMax; }
    else if (mDetectionPropertiesNew.v.thresholdCircumferenceChange < mDetectionPropertiesNew.p.circumferenceChangeThreshold)
    {        mDetectionPropertiesNew.v.thresholdCircumferenceChange = mDetectionPropertiesNew.p.circumferenceChangeThreshold; }

    if (mDetectionProperties.v.thresholdDisplacementChange > imageWdth)
    {   mDetectionProperties.v.thresholdDisplacementChange = imageWdth; }
    else if (mDetectionProperties.v.thresholdDisplacementChange < mDetectionProperties.p.displacementChangeThreshold)
    {        mDetectionProperties.v.thresholdDisplacementChange = mDetectionProperties.p.displacementChangeThreshold; }

    if (mDetectionPropertiesNew.v.curvatureOffset > 180)
    {   mDetectionPropertiesNew.v.curvatureOffset = 180; }
    else if (mDetectionPropertiesNew.v.curvatureOffset < mDetectionPropertiesNew.p.curvatureOffset)
    {        mDetectionPropertiesNew.v.curvatureOffset = mDetectionPropertiesNew.p.curvatureOffset; }

    if (mDetectionPropertiesNew.v.priorCertainty > certaintyUpperLimit)
    {   mDetectionPropertiesNew.v.priorCertainty = certaintyUpperLimit; }
    else if (mDetectionPropertiesNew.v.priorCertainty < certaintyLowerLimit)
    {        mDetectionPropertiesNew.v.priorCertainty = certaintyLowerLimit; }

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
