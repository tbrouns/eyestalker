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
    int gradientWindowLength = glintAOI.wdth;
    int glintRadius = round(0.5 * glintAOI.wdth);

    uchar *ptr = img.data;
    
    std::vector<double> imageGradient(searchAOI.wdth * searchAOI.hght, 0.0);
    
    std::vector<int> dZ(8);
    dZ[0] = -1;
    dZ[1] = -searchAOI.wdth - 1;
    dZ[2] = -searchAOI.wdth;
    dZ[3] = -searchAOI.wdth + 1;
    dZ[4] = 1;
    dZ[5] = searchAOI.wdth + 1;
    dZ[6] = searchAOI.wdth;
    dZ[7] = searchAOI.wdth - 1;
    
    for (int y = gradientWindowLength; y < searchAOI.hght - gradientWindowLength; y++)
    {
        for (int x = gradientWindowLength; x < searchAOI.wdth - gradientWindowLength; x++)
        {
            int i = searchAOI.wdth * y + x; // gradient coordinates
            int j = imgWidth * (y + searchAOI.yPos) + (x + searchAOI.xPos); // image coordinates
            
            double centreSum = ptr[j];
            for (int m = 0; m < 8; m++) { centreSum += ptr[j + dZ[m]]; }
            
            double surroundSum = 0;
            for (int m = 0; m < 8; m++) { surroundSum += ptr[j + gradientWindowLength * dZ[m]]; }
            imageGradient[i] = centreSum / surroundSum;
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
    innerAOI.xPos = 0;
    innerAOI.yPos = 0;
    
    double minPupilIntensity = pow(10, 10); // arbitrarily large value
    
    int innerAOIArea = (innerAOI.wdth - 1) * (innerAOI.hght - 1);
    
    for (int y = 0; y < searchAOI.hght - innerAOI.hght; y++)
    {
        for (int x = 0; x < searchAOI.wdth - innerAOI.wdth; x++)
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
                int glintTopLeftIndex  = searchAOI.wdth * glintAOI.yPos + glintAOI.xPos;
                int glintTopRghtIndex  = searchAOI.wdth * glintAOI.yPos + glintAOI.xPos + glintAOI.wdth;
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

std::vector<int> radialGradient(const cv::Mat& img, int kernelSize, double pupilXCentre, double pupilYCentre)
{
    int numDirections = 8;

    std::vector<int> dX(numDirections);
    dX[0] =  1;
    dX[1] =  1;
    dX[2] =  0;
    dX[3] = -1;
    dX[4] = -1;
    dX[5] = -1;
    dX[6] =  0;
    dX[7] =  1;

    std::vector<int> dY(numDirections);
    dY[0] =  0;
    dY[1] = -1;
    dY[2] = -1;
    dY[3] = -1;
    dY[4] =  0;
    dY[5] =  1;
    dY[6] =  1;
    dY[7] =  1;

    uchar *ptr = img.data;
    int width  = img.cols;
    int height = img.rows;

    std::vector<int> gradientVector(width * height, 0);

    int kernelRadius = (kernelSize - 1) / 2;

    for (int y = kernelRadius; y < height - kernelRadius; y++)
    {
        for (int x = kernelRadius; x < width - kernelRadius; x++)
        {
            double dx = x - pupilXCentre;
            double dy = pupilYCentre - y;

            double theta;
            if (dx != 0 && dy != 0)
            {
                theta = atan2(dy,dx);
                if (theta < 0) { theta = theta + 2 * M_PI; }
            }
            else if (dx == 0 && dy != 0)
            {
                if (dy > 0) { theta = 0.5 * M_PI; }
                if (dy < 0) { theta = 1.5 * M_PI; }
            }
            else if (dx != 0 && dy == 0)
            {
                if (dx > 0) { theta = 0; }
                if (dx < 0) { theta = M_PI; }
            }
            else { theta = 0; } // arbitrary choice

            int positiveDirection = round(theta * (numDirections - 1) / (2 * M_PI)); // direction index of gradient
            int negativeDirection = (int)(positiveDirection + 0.5 * numDirections) % numDirections;

            int neighbourIndexPositive = width * (y + dY[positiveDirection] * kernelRadius) + (x + dX[positiveDirection] * kernelRadius);
            int neighbourIndexNegative = width * (y + dY[negativeDirection] * kernelRadius) + (x + dX[negativeDirection] * kernelRadius);

            double val = ptr[neighbourIndexPositive] - ptr[neighbourIndexNegative];
            if (val < 0) { val = 0; }
            int centreIndex = width * y + x;
            gradientVector[centreIndex] = val;
        }
    }

    return gradientVector;
}

std::vector<int> nonMaximumSuppresion(const std::vector<int>& gradient, int width, int height, double pupilXCentre, double pupilYCentre)
{
    std::vector<int> gradientSuppressed = gradient;

    int kernelPerimeter = 8;

    std::vector<int> dX(kernelPerimeter);
    dX[0] =  1;
    dX[1] =  1;
    dX[2] =  0;
    dX[3] = -1;
    dX[4] = -1;
    dX[5] = -1;
    dX[6] =  0;
    dX[7] =  1;

    std::vector<int> dY(kernelPerimeter);
    dY[0] =  0;
    dY[1] = -1;
    dY[2] = -1;
    dY[3] = -1;
    dY[4] =  0;
    dY[5] =  1;
    dY[6] =  1;
    dY[7] =  1;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int iPixel = y * width + x;
            if (gradient[iPixel] == 0) { continue; }

            double dx = x - pupilXCentre;
            double dy = pupilYCentre - y;

            double theta;
            if (dx != 0 && dy != 0)
            {
                theta = atan2(dy,dx);
                if (theta < 0) { theta = theta + 2 * M_PI; }
            }
            else { theta = 0; }   // arbitrary choice

            int i = round(theta * kernelPerimeter / (2 * M_PI));
            int j = (kernelPerimeter / 2) + i;

            i = i % kernelPerimeter;
            j = j % kernelPerimeter;

            int iNeighbour = width * (y + dY[i]) + (x + dX[i]);
            int jNeighbour = width * (y + dY[j]) + (x + dX[j]);

            if (gradient[iPixel] < gradient[iNeighbour] || gradient[iPixel] < gradient[jNeighbour])
            {
                gradientSuppressed[iPixel] = 0;
            }
        }
    }

    return gradientSuppressed;
}

std::vector<int> hysteresisTracking(const std::vector<int>& gradientSuppressed, int width, int height, int thresholdHigh, int thresholdLow)
{
    std::vector<int> dX(8);
    dX[0] = -1;
    dX[1] = -1;
    dX[2] =  0;
    dX[3] =  1;
    dX[4] =  1;
    dX[5] =  1;
    dX[6] =  0;
    dX[7] = -1;

    std::vector<int> dY(8);
    dY[0] =  0;
    dY[1] = -1;
    dY[2] = -1;
    dY[3] = -1;
    dY[4] =  0;
    dY[5] =  1;
    dY[6] =  1;
    dY[7] =  1;

    int imgSize = width * height;
    std::vector<int> edgePointVector(imgSize, 0);

    for (int iPixel = 0; iPixel < imgSize; iPixel++)
    {
        if (gradientSuppressed[iPixel] >= thresholdHigh && edgePointVector[iPixel] == 0)
        {
            edgePointVector[iPixel] = 1;

            std::vector<int> oldEdgeIndices(1);
            oldEdgeIndices[0] = iPixel;

            do
            {
                std::vector<int> newEdgeIndices;
                int numIndices = oldEdgeIndices.size();

                for (int iEdgePoint = 0; iEdgePoint < numIndices; iEdgePoint++)
                {
                    int centreIndex = oldEdgeIndices[iEdgePoint];
                    int centreXPos  = centreIndex % width;
                    int centreYPos  = (centreIndex - centreXPos) / width;

                    for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
                    {
                        int neighbourXPos = centreXPos + dX[m];
                        int neighbourYPos = centreYPos + dY[m];

                        if (neighbourXPos < 0 || neighbourXPos >= width ||
                                neighbourYPos < 0 || neighbourYPos >= height)
                        { continue; } // neighbour is out-of-bounds

                        int neighbourIndex = width * neighbourYPos + neighbourXPos;

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

std::vector<int> sharpenEdges(std::vector<int>& binaryImageVectorRaw, int haarWidth, int haarHeight)
{
    std::vector<int> binaryImageVector = binaryImageVectorRaw;

    std::vector<int> dX(8);
    dX[0] =  0;
    dX[1] =  1;
    dX[2] = -1;
    dX[3] =  1;
    dX[4] =  0;
    dX[5] = -1;
    dX[6] =  1;
    dX[7] = -1;

    std::vector<int> dY(8);
    dY[0] = -1;
    dY[1] =  1;
    dY[2] =  0;
    dY[3] = -1;
    dY[4] =  1;
    dY[5] = -1;
    dY[6] =  0;
    dY[7] =  1;

    for (int yCentre = 0; yCentre < haarHeight; yCentre++)
    {
        for (int xCentre = 0; xCentre < haarWidth; xCentre++)
        {
            int iCentre = haarWidth * yCentre + xCentre;

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

                        if (xNeighbour < 0 || xNeighbour >= haarWidth ||yNeighbour < 0 || yNeighbour >= haarHeight)
                        {
                            continue; // neighbour is out-of-bounds
                        }

                        int iNeighbour = haarWidth * yNeighbour + xNeighbour;

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
                        int iOpposite = haarWidth * yOpposite + xOpposite;

                        if
                                (xOpposite < 0 || xOpposite >= haarWidth  || // ... AND opposite pixel is out-of-bounds
                                 yOpposite < 0 || yOpposite >= haarHeight ||
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

std::vector<int> cannyConversion(const cv::Mat& img, int haarWidth, int haarHeight)
{
    int haarSize = haarWidth * haarHeight;
    uchar *ptr_img = img.data;
    std::vector<int> binaryImageVectorRaw(haarSize);

    for (int i = 0; i < haarSize; i++)
    {
        if (ptr_img[i] == 255)  { binaryImageVectorRaw[i] = 1; }
        else                    { binaryImageVectorRaw[i] = 0; }
    }

    return binaryImageVectorRaw;
}

std::vector<int> getEdgeIndices(const std::vector<int>& binaryImageVector)
{
    int haarSize = binaryImageVector.size();
    std::vector<int> cannyEdgeIndices;
    for (int i = 0; i < haarSize; i++)
    { if (binaryImageVector[i] == 1) { cannyEdgeIndices.push_back(i); }}
    return cannyEdgeIndices;
}

std::vector<edgeProperties> edgeSegmentation(const cv::Mat& img, std::vector<int>& p, int haarWidth, int haarHeight, double pupilXCentre, double pupilYCentre, double curvatureLowerLimit, double curvatureUpperLimit, bool USE_PRIOR_INFORMATION)
{
    std::vector<edgeProperties> vEdgePropertiesAll; // new structure containing length and indices of all edges

    // scanned neighbours
    std::vector<int> dX(8);
    dX[0] = -1;
    dX[1] = -1;
    dX[2] =  0;
    dX[3] =  1;
    dX[4] =  1;
    dX[5] =  1;
    dX[6] =  0;
    dX[7] = -1;

    std::vector<int> dY(8);
    dY[0] =  0;
    dY[1] = -1;
    dY[2] = -1;
    dY[3] = -1;
    dY[4] =  0;
    dY[5] =  1;
    dY[6] =  1;
    dY[7] =  1;

    std::vector<int> dZ(8);
    dZ[0] = -1;
    dZ[1] = -haarWidth - 1;
    dZ[2] = -haarWidth;
    dZ[3] = -haarWidth + 1;
    dZ[4] =  1;
    dZ[5] =  haarWidth + 1;
    dZ[6] =  haarWidth;
    dZ[7] =  haarWidth - 1;
    
    int startEdgeIndex = 0;

    while (1)
    {
        bool NEW_EDGE_FOUND = false;

        if (USE_PRIOR_INFORMATION)
        {
            for (int m = 0; m < 8 && !NEW_EDGE_FOUND; m++) // Starburst
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

                        if (x < 0 || x >= haarWidth || y < 0 || y >= haarHeight) { STOP_SEARCH = true; break; }

                        int i = y * haarWidth + x;

                        if (p[i] == 1)
                        {
                            startEdgeIndex = i;
                            NEW_EDGE_FOUND = true;
                            STOP_SEARCH    = true;
                        }
                        else if (p[i] > 1) { STOP_SEARCH = true; }
                    }
                }
            }
        }
        else
        {
            for (int i = startEdgeIndex; i < haarWidth * haarHeight && !NEW_EDGE_FOUND; i++)
            {
                if (p[i] == 1)
                {
                    startEdgeIndex = i;
                    NEW_EDGE_FOUND = true;
                }
            }
        }

        if (!NEW_EDGE_FOUND) { break; } // no (more) edges found

        // find edge endpoint
        
        std::vector<int> allEdgeIndicesRaw(1); // all indices belonging to the edge, out-of-order
        allEdgeIndicesRaw[0] = startEdgeIndex;
        
        p[startEdgeIndex] = 2; // tag pixel
        
        {
            std::vector<int> rawEdgeIndices(1);
            rawEdgeIndices[0] = startEdgeIndex;
            
            bool edgePointFound = true;
            
            while (edgePointFound)
            {
                std::vector<int> newEdgePoints;
                
                edgePointFound = false;
                
                for (int iEdgePoint = 0, numberOfNewPoints = rawEdgeIndices.size(); iEdgePoint < numberOfNewPoints; iEdgePoint++) // loop through all newly added unchecked edge points
                {
                    int centreIndex = rawEdgeIndices[iEdgePoint]; // index of current edge point

                    int centreXPos = centreIndex % haarWidth;
                    int centreYPos = (centreIndex - centreXPos) / haarWidth;

                    for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
                    {
                        int neighbourXPos = centreXPos + dX[m];
                        int neighbourYPos = centreYPos + dY[m];
                        
                        if (neighbourXPos < 0 || neighbourXPos >= haarWidth || neighbourYPos < 0 || neighbourYPos >= haarHeight)
                        {
                            continue; // neighbour is out-of-bounds
                        }
                        
                        int neighbourIndex = haarWidth * neighbourYPos + neighbourXPos;

                        if (p[neighbourIndex] == 1) // if neighbouring point is filled ...
                        {
                            p[neighbourIndex] = 2; // ... then tag it
                            
                            newEdgePoints.push_back(neighbourIndex); // edge points to-be-checked
                            allEdgeIndicesRaw.push_back(neighbourIndex); // all edge points
                            
                            edgePointFound = true;
                        }
                    }
                }
                
                rawEdgeIndices = newEdgePoints;
                newEdgePoints.clear();
            }
        }

        int edgeStartIndex = allEdgeIndicesRaw.back(); // edge terminal
        p[edgeStartIndex] = 3; // tag edge terminal
        
        std::vector<int> allEdgeIndices; // all edge indices, in sequence
        
        std::vector<int> branchStartIndices(1);
        branchStartIndices[0] = edgeStartIndex;
        
        int numbranches = 1; // start off with one branch

        while (numbranches > 0)
        {
            std::vector<int> branchSizes;                            // record length of each branch
            std::vector<std::vector<int>> branchIndicesVectors;      // record indices of each branch
            std::vector<std::vector<int>> branchStartIndicesVectors; // record start index of each branch

            for (int iBranch = 0; iBranch < numbranches; iBranch++) // stops when junction is encountered
            {
                int branchStartIndex = branchStartIndices[iBranch];
                int currentIndex = branchStartIndex;

                std::vector<int> branchIndices(1); // record all indices of branch
                branchIndices[0] = branchStartIndex;
                
                std::vector<int> branchStartIndicesNew;

                while(1)
                {
                    std::vector<int> newBranchIndices;
                    
                    int centreIndex = currentIndex;
                    int centreXPos = centreIndex % haarWidth;
                    int centreYPos = (centreIndex - centreXPos) / haarWidth;

                    for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
                    {
                        int neighbourXPos = centreXPos + dX[m];
                        int neighbourYPos = centreYPos + dY[m];

                        if (neighbourXPos < 0 || neighbourXPos >= haarWidth || neighbourYPos < 0 || neighbourYPos >= haarHeight)
                        {
                            continue; // neighbour is out-of-bounds
                        }

                        int neighbourIndex = haarWidth * neighbourYPos + neighbourXPos;
                        
                        if (p[neighbourIndex] == 2) // if neighbouring point was tagged previously ...
                        {
                            p[neighbourIndex] = 3; // ... give it a new tag
                            
                            newBranchIndices.push_back(neighbourIndex);
                            currentIndex = neighbourIndex;
                        }
                    }
                    
                    int numberOfNewPoints = newBranchIndices.size();

                    if (numberOfNewPoints == 1) // only one new point is added: continuation of single branch
                    {
                        branchIndices.push_back(newBranchIndices[0]);
                    }
                    else if (numberOfNewPoints >= 2) // multiple points are added: edge branches off
                    {
                        branchStartIndicesNew = newBranchIndices;
                        break;
                    }
                    else // no new points added
                    {
                        break;
                    }
                }

                branchSizes.push_back(branchIndices.size());
                branchIndicesVectors.push_back(branchIndices);
                branchStartIndicesVectors.push_back(branchStartIndicesNew);
            }

            int maxIndex = std::distance(branchSizes.begin(), std::max_element(branchSizes.begin(), branchSizes.end()));
            allEdgeIndices.insert(std::end(allEdgeIndices), std::begin(branchIndicesVectors[maxIndex]), std::end(branchIndicesVectors[maxIndex])); // only add longest branch

            branchStartIndices = branchStartIndicesVectors[maxIndex];
            numbranches = branchStartIndices.size();
        }

        // give pixels that have been added to new outline a new tag

        for (int iEdgePoint = 0, edgeLength = allEdgeIndices.size(); iEdgePoint < edgeLength; iEdgePoint++)
        {
            p[allEdgeIndices[iEdgePoint]] = 4;
        }

        // remove tag from pixels that have been tagged before, but not included in new outline
        
        for (int iEdgePoint = 0, edgeLength = allEdgeIndicesRaw.size(); iEdgePoint < edgeLength; iEdgePoint++)
        {
            int edgePointIndex = allEdgeIndicesRaw[iEdgePoint];
            
            if (p[edgePointIndex] == 2 || p[edgePointIndex] == 3)
            {
                p[edgePointIndex] = 1;
            }
        }

        int edgeLength = allEdgeIndices.size();

        if (edgeLength < curvatureWindowLength) // ignore short edges
        {
            for (int iEdgePoint = 0; iEdgePoint < edgeLength; iEdgePoint++)
            {
                int edgePointIndex = allEdgeIndices[iEdgePoint];
                p[edgePointIndex]  = -1; // make them transparent to starburst
            }

            continue;
        }

        // calculate curvature

        // calculate directions of edge points
        
        std::vector<double> xOrientation(8);
        std::vector<double> yOrientation(8);
        
        xOrientation[0] = -1.0;
        yOrientation[0] =  0.0;
        xOrientation[1] = -sqrt(0.5);
        yOrientation[1] = -sqrt(0.5);
        xOrientation[2] =  0.0;
        yOrientation[2] = -1.0;
        xOrientation[3] =  sqrt(0.5);
        yOrientation[3] = -sqrt(0.5);
        xOrientation[4] = 1.0;
        yOrientation[4] = 0.0;
        xOrientation[5] = sqrt(0.5);
        yOrientation[5] = sqrt(0.5);
        xOrientation[6] = 0.0;
        yOrientation[6] = 1.0;
        xOrientation[7] = -sqrt(0.5);
        yOrientation[7] =  sqrt(0.5);
        
        std::vector<double> edgeXOrientations(edgeLength); // direction vector for each edge pixel
        std::vector<double> edgeYOrientations(edgeLength);
        
        for (int iEdgePoint = 0; iEdgePoint < edgeLength - 1; iEdgePoint++)
        {
            int centreIndex = allEdgeIndices[iEdgePoint];
            int nextCentreIndex = allEdgeIndices[iEdgePoint + 1];
            
            for (int m = 0; m < 8; m++)
            {
                int adjacentIndex = centreIndex + dZ[m];
                
                if (nextCentreIndex == adjacentIndex)
                {
                    edgeXOrientations[iEdgePoint] = xOrientation[m];
                    edgeYOrientations[iEdgePoint] = yOrientation[m];
                    break;
                }
            }
        }
        
        // add last indices and give them directions of second-to-last
        
        edgeXOrientations[edgeLength - 1] = edgeXOrientations[edgeLength - 2];
        edgeYOrientations[edgeLength - 1] = edgeYOrientations[edgeLength - 2];

        // calculate curvature and normal lines to curve

        std::vector<double> edgePointCurvatures(edgeLength, 360.0);
        std::vector<double> edgePointXNormals(edgeLength, 0.0);
        std::vector<double> edgePointYNormals(edgeLength, 0.0);

        for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeLength - curvatureWindowLength; iEdgePoint++)
        {
            // calculate tangents

            // first window

            double firstXTangent = calculateMean(std::vector<double>(&edgeXOrientations[iEdgePoint - curvatureWindowLength],&edgeXOrientations[iEdgePoint - 1]));
            double firstYTangent = calculateMean(std::vector<double>(&edgeYOrientations[iEdgePoint - curvatureWindowLength],&edgeYOrientations[iEdgePoint - 1]));

            // second window

            double secndXTangent = calculateMean(std::vector<double>(&edgeXOrientations[iEdgePoint + 1],&edgeXOrientations[iEdgePoint + curvatureWindowLength]));
            double secndYTangent = calculateMean(std::vector<double>(&edgeYOrientations[iEdgePoint + 1],&edgeYOrientations[iEdgePoint + curvatureWindowLength]));

            // calculate vector difference

            double vectorAngle = atan2(secndYTangent, secndXTangent) - atan2(firstYTangent, firstXTangent);

            if      (vectorAngle >  M_PI) { vectorAngle = vectorAngle - 2 * M_PI; }
            else if (vectorAngle < -M_PI) { vectorAngle = vectorAngle + 2 * M_PI; }

            edgePointCurvatures[iEdgePoint] = 180 * vectorAngle / M_PI; // in degrees

            edgePointXNormals[iEdgePoint] = -firstXTangent + secndXTangent;
            edgePointYNormals[iEdgePoint] = -firstYTangent + secndYTangent;
        }

        // find majority sign of curvature
        
        int edgeCurvatureSign = 1;
        
        {
            int numberOfPositiveCurvatures = 0;
            int numberOfNegativeCurvatures = 0;
            
            for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeLength - curvatureWindowLength; iEdgePoint++)
            {
                double   edgePointCurvature = edgePointCurvatures[iEdgePoint];
                if      (edgePointCurvature > 0) { numberOfPositiveCurvatures++; }
                else if (edgePointCurvature < 0) { numberOfNegativeCurvatures++; }
            }
            
            if (numberOfPositiveCurvatures >= numberOfNegativeCurvatures) { edgeCurvatureSign =  1; }
            else                                                          { edgeCurvatureSign = -1; }
        }

        // find breakpoints
        
        std::vector<int> breakPointIndices; // position of breakpoints
        breakPointIndices.push_back(0); // add first point

        for (int iEdgePoint = curvatureWindowLength; iEdgePoint < edgeLength - curvatureWindowLength; iEdgePoint++)
        {
            double edgePointCurvature = edgePointCurvatures[iEdgePoint];

            if ((std::abs(edgePointCurvature) >= curvatureUpperLimit) || (edgeCurvatureSign * edgePointCurvature <= curvatureLowerLimit))
            {
                breakPointIndices.push_back(iEdgePoint);
            }
        }

        breakPointIndices.push_back(edgeLength - 1); // add last point
        
        // evaluate each partial edge
        
        std::sort(breakPointIndices.begin(), breakPointIndices.end());
        
        int numberOfBreakPoints = breakPointIndices.size();
        
        for (int iBreakPoint = 0; iBreakPoint < numberOfBreakPoints - 1; iBreakPoint++)
        {
            int iStartBreakPoint  = breakPointIndices[iBreakPoint] + 1;
            int iEndBreakPoint    = breakPointIndices[iBreakPoint + 1];
            int partialEdgeLength = iEndBreakPoint - iStartBreakPoint;

            // grab indices of edge points

            if (partialEdgeLength < curvatureWindowLength) { continue; } // ignore short edges
            std::vector<int> partialEdgeIndices(allEdgeIndices.begin() + iStartBreakPoint, allEdgeIndices.begin() + iEndBreakPoint);

            // calculate pixel intensities within inner curve of edge points

            std::vector<double> partialEdgeIntensities(partialEdgeLength);

            {
                uchar *ptr_img = img.data;

                for (int iEdgePoint = 0; iEdgePoint < partialEdgeLength; iEdgePoint++)
                {
                    int edgePointIndex = partialEdgeIndices[iEdgePoint];
                    int edgePointXPos  = edgePointIndex % haarWidth;
                    int edgePointYPos  = (edgePointIndex - edgePointXPos) / haarWidth;

                    int offsetXPos = edgePointXPos + edgeIntensityPositionOffset * ceil2(edgePointXNormals[iEdgePoint]);
                    int offsetYPos = edgePointYPos + edgeIntensityPositionOffset * ceil2(edgePointYNormals[iEdgePoint]);

                    if (offsetXPos < 0 || offsetXPos > haarWidth || offsetYPos < 0 || offsetYPos > haarHeight)
                    {       partialEdgeIntensities[iEdgePoint] = ptr_img[edgePointXPos + edgePointYPos * haarWidth]; }
                    else {  partialEdgeIntensities[iEdgePoint] = ptr_img[   offsetXPos +    offsetYPos * haarWidth]; }
                }
            }

            double avgIntensity = calculateMean(partialEdgeIntensities); // calculate average intensity

            // calculate average absolute curvature

            double avgCurvature = 0;
            double maxCurvature;
            double minCurvature;

            {
                std::vector<double> partialEdgeCurvatures(edgePointCurvatures.begin() + iStartBreakPoint, edgePointCurvatures.begin() + iEndBreakPoint);
                std::vector<double> partialEdgeCurvaturesNew;

                for (int iEdgePoint = 0; iEdgePoint < partialEdgeLength; iEdgePoint++)
                {
                    double curvature = partialEdgeCurvatures[iEdgePoint];
                    if (curvature < 180.0) { partialEdgeCurvaturesNew.push_back(curvature); }
                }

                int partialEdgeLengthNew = partialEdgeCurvaturesNew.size();

                if (partialEdgeLengthNew > 0)
                {
                    for (int iEdgePoint = 0; iEdgePoint < partialEdgeLengthNew; iEdgePoint++)
                    {
                        avgCurvature += std::abs(partialEdgeCurvaturesNew[iEdgePoint] / partialEdgeLengthNew);
                    }

                    maxCurvature = *std::max_element(std::begin(partialEdgeCurvaturesNew), std::end(partialEdgeCurvaturesNew));
                    minCurvature = *std::min_element(std::begin(partialEdgeCurvaturesNew), std::end(partialEdgeCurvaturesNew));
                }
                else
                {
                    avgCurvature = 360;
                    maxCurvature = 360;
                    minCurvature = 360;
                }
            }

            // add additional adjacent indices that were removed by morphological operation
            
            for (int iEdgePoint = 0; iEdgePoint < partialEdgeLength; iEdgePoint++)
            {
                int centreIndex = partialEdgeIndices[iEdgePoint];

                int centreXPos = centreIndex % haarWidth;
                int centreYPos = (centreIndex - centreXPos) / haarWidth;

                for (int m = 0; m < 8; m++) // loop through 8-connected environment
                {
                    int neighbourXPos = centreXPos + dX[m];
                    int neighbourYPos = centreYPos + dY[m];

                    if (neighbourXPos < 0 || neighbourXPos >= haarWidth || neighbourYPos < 0 || neighbourYPos >= haarHeight)
                    {
                        continue; // neighbour is out-of-bounds
                    }

                    int neighbourIndex = haarWidth * neighbourYPos + neighbourXPos;

                    if (p[neighbourIndex] == -1) // if neighbouring point was canny edge point that was removed by morphological operation then ...
                    {
                        p[neighbourIndex] = 4; // ... tag it and ...
                        partialEdgeIndices.push_back(neighbourIndex); // ... add it to the (partial) edge
                    }
                }
            }

            int numEdgePoints = partialEdgeIndices.size();

            edgeProperties mEdgeProperties;
            mEdgeProperties.curvatureAvg = avgCurvature;
            mEdgeProperties.curvatureMax = maxCurvature;
            mEdgeProperties.curvatureMin = minCurvature;
            mEdgeProperties.intensity    = avgIntensity;
            mEdgeProperties.length       = partialEdgeLength;
            mEdgeProperties.size         = numEdgePoints;
            mEdgeProperties.pointIndices = partialEdgeIndices;

            vEdgePropertiesAll.push_back(mEdgeProperties);
        }
    }

    return vEdgePropertiesAll; // return structure
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

    //#ifdef __linux__
    //        #pragma omp parallel for
    //#endif
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
                    combiEdgeIndices[jEdge]      = vEdgePropertiesAll[iEdge].index;
                    combiEdgeIntensities[jEdge]  = vEdgePropertiesAll[iEdge].intensity;
                    combiEdgeLengths[jEdge]      = vEdgePropertiesAll[iEdge].length;
                    combiEdgeSizes[jEdge]        = vEdgePropertiesAll[iEdge].size;
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
        cv::Mat imagePupilBGR = imageOriginalBGR(outerRect);

        // Convert back to grayscale

        cv::Mat imagePupilGray;
        cv::cvtColor(imagePupilBGR, imagePupilGray, cv::COLOR_BGR2GRAY);

        ///////////////////////////////////////////////////////////////////////
        /////////////////////// CANNY EDGE DETECTION  /////////////////////////
        ///////////////////////////////////////////////////////////////////////       

        cv::Mat imagePupilGrayBlurred;
        int cannyBlurLevel = 2 * mDetectionProperties.p.cannyBlurLevel - 1; // should be odd
        if (cannyBlurLevel > 0) { cv::GaussianBlur(imagePupilGray, imagePupilGrayBlurred, cv::Size(cannyBlurLevel, cannyBlurLevel), 0, 0);
        } else                  { imagePupilGrayBlurred = imagePupilGray; }

        double xPosPredictionRelative = mDetectionProperties.v.xPosPrediction - outerAOI.xPos;
        double yPosPredictionRelative = mDetectionProperties.v.yPosPrediction - outerAOI.yPos;

        std::vector<int> cannyEdges; // binary vector

        if (USE_PRIOR_INFORMATION)
        {
            std::vector<int>  imgGradient           = radialGradient(imagePupilGrayBlurred, mDetectionProperties.p.cannyKernelSize, xPosPredictionRelative, yPosPredictionRelative);
            std::vector<int>  imgGradientSuppressed = nonMaximumSuppresion(imgGradient, outerAOI.wdth, outerAOI.hght, xPosPredictionRelative, yPosPredictionRelative);
            cannyEdges = hysteresisTracking(imgGradientSuppressed, outerAOI.wdth, outerAOI.hght, mDetectionProperties.p.cannyThresholdHigh, mDetectionProperties.p.cannyThresholdLow);
        }
        else
        {
            cv::Mat imageCannyEdges;
            cv::Canny(imagePupilGrayBlurred, imageCannyEdges, 30 * mDetectionProperties.p.cannyThresholdHigh, 30 * mDetectionProperties.p.cannyThresholdLow, 5);
            cannyEdges = cannyConversion(imageCannyEdges, outerAOI.wdth, outerAOI.hght);
        }

        std::vector<int> edgeIndices         = getEdgeIndices(cannyEdges);
        std::vector<int> cannyEdgesSharpened = sharpenEdges(cannyEdges, outerAOI.wdth, outerAOI.hght); // Morphological operation

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

        std::vector<edgeProperties> vEdgePropertiesAll = edgeSegmentation(imagePupilGray, cannyEdgesSharpened, outerAOI.wdth, outerAOI.hght, xPosPredictionRelative, yPosPredictionRelative, curvatureLowerLimit, curvatureUpperLimit, USE_PRIOR_INFORMATION);

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

            mEllipseProperties = findBestEllipseFit(vEdgePropertiesNew, outerAOI, mDetectionProperties, USE_PRIOR_INFORMATION); // ellipse fitting

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
            cv::Mat imagePupilBGR = imageOriginalBGR(pupilROI);
            cv::Mat imagePupilGray;
            cv::cvtColor(imagePupilBGR, imagePupilGray, cv::COLOR_BGR2GRAY);

            mDetectionPropertiesNew.m.imagePupil = imagePupilGray;
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
