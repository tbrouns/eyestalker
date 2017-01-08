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

#include "headers/eyetracking.h"

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

std::vector<unsigned int> calculateIntImg(const cv::Mat& img, int imgWidth, int startX, int startY, int width, int height)
{
    uchar *ptr = img.data;
    
    std::vector<unsigned int> integralImage(width * height); // unsigned due to large positive values
    
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int i = width * y + x; // integral image coordinates
            int j = imgWidth * (y + startY) + (x + startX); // image coordinates
            
            double val = ptr[j];

            if (x == 0 && y == 0) { integralImage[i] = val; } // first point
            else if (y == 0)      { integralImage[i] = val + integralImage[i - 1]; } // first row
            else if (x == 0)      { integralImage[i] = val + integralImage[i - width]; } // first column
            else                  { integralImage[i] = val + integralImage[i - 1] + integralImage[i - width] - integralImage[i - width - 1]; }
        }
    }
    
    return integralImage;
}

haarProperties detectGlint(const cv::Mat& img, int imgWidth, int startX, int startY, int width, int height, int glintSize)
{
    int gradientWindowLength = glintSize;
    int glintRadius = round(0.5 * glintSize);

    uchar *ptr = img.data;
    
    std::vector<double> imageGradient(width * height, 0.0);
    
    std::vector<int> dZ(8);
    dZ[0] = -1;
    dZ[1] = -width - 1;
    dZ[2] = -width;
    dZ[3] = -width + 1;
    dZ[4] = 1;
    dZ[5] = width + 1;
    dZ[6] = width;
    dZ[7] = width - 1;
    
    for (int y = gradientWindowLength; y < height - gradientWindowLength; y++)
    {
        for (int x = gradientWindowLength; x < width - gradientWindowLength; x++)
        {
            int i = width * y + x; // gradient coordinates
            int j = imgWidth * (y + startY) + (x + startX); // image coordinates
            
            double centre = ptr[j];
            
            for (int m = 0; m < 8; m++)
            {
                centre += ptr[j + dZ[m]];
            }
            
            double surround = 0;
            
            for (int m = 0; m < 8; m++)
            {
                surround += ptr[j + gradientWindowLength * dZ[m]];
            }
            
            imageGradient[i] = centre / surround;
        }
    }
    
    int glintIndex = std::distance(imageGradient.begin(), std::max_element(imageGradient.begin(), imageGradient.end()));
    
    int x = glintIndex % width;
    int y = (glintIndex - x) / width;
    
    haarProperties glint;
    
    glint.xPos = x - glintRadius;
    glint.yPos = y - glintRadius;
    
    return glint;
}

haarProperties detectPupilApprox(const std::vector<unsigned int>& I, int width, int height, int haarWidth, int haarHeight, int glintX, int glintY, int glintSize)
{
    int pupilX = 0;
    int pupilY = 0;
    
    double minPupilIntensity = pow(10, 10); // arbitrarily large value
    
    int pupilArea = (haarWidth - 1) * (haarHeight - 1);
    
    for (int y = 0; y < height - haarHeight; y++)
    {
        for (int x = 0; x < width - haarWidth; x++)
        {
            // vertices of inner square
            
            int topLeftX = x;
            int topLeftY = y;
            
            int backRightX = topLeftX + (haarWidth - 1);
            int backRightY = topLeftY + (haarHeight - 1);
            
            int topLeftIndex  = width * topLeftY + topLeftX;
            int topRghtIndex  = topLeftIndex + (haarWidth - 1);
            int backLeftIndex = topLeftIndex + (haarHeight - 1) * width;
            int backRghtIndex = topRghtIndex + backLeftIndex - topLeftIndex;
            
            // calculate glint intensity
            
            double glintIntensity = 0.0;
            double glintArea = 0.0;
            
            bool glintWithinHaarDetector = false; // flag for glint overlap
            
            std::vector<int> z(4);
            z[0] = 0;
            z[1] = 0;
            z[2] = glintSize;
            z[3] = glintSize;
            
            // check if glint overlaps with Haar detector
            
            for (int m = 0; m < 4; m++)
            {
                int n = (m + 1) % 4;
                
                if (glintX + z[m] >= topLeftX && glintX + z[m] <= backRightX)
                {
                    if (glintY + z[n] >= topLeftY && glintY + z[n] <= backRightY)
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
                
                if (glintX < topLeftX)
                {
                    glintOverlapsLeftEdge = true;
                }
                else if (glintX + glintSize > backRightX)
                {
                    glintOverlapsRightEdge = true;
                }
                
                if (glintY < topLeftY)
                {
                    glintOverlapsTopEdge = true;
                }
                else if (glintY + glintSize > backRightY)
                {
                    glintOverlapsBottomEdge = true;
                }
                
                // coordinates of corners of glint square
                int glintTopLeftIndex  = width * glintY + glintX;
                int glintTopRghtIndex  = width * glintY + glintX + glintSize;
                int glintBackLeftIndex = width * (glintY + glintSize) + glintX;
                int glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                
                // check if glint square overlaps with edge or corner of pupil square
                
                // check edge overlap
                
                if (!glintOverlapsLeftEdge && !glintOverlapsRightEdge && glintOverlapsTopEdge) // top edge
                {
                    glintTopLeftIndex  = width * topLeftY + glintX;
                    glintTopRghtIndex  = width * topLeftY + glintX + glintSize;
                    glintBackLeftIndex = width * (glintY + glintSize) + glintX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                
                if (!glintOverlapsLeftEdge && !glintOverlapsRightEdge && glintOverlapsBottomEdge) // bottom edge
                {
                    glintTopLeftIndex  = width * glintY + glintX;
                    glintTopRghtIndex  = width * glintY + glintX + glintSize;
                    glintBackLeftIndex = width * backRightY + glintX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                
                if (glintOverlapsLeftEdge && !glintOverlapsTopEdge && !glintOverlapsBottomEdge) // left edge
                {
                    glintTopLeftIndex  = width * glintY + topLeftX;
                    glintTopRghtIndex  = width * glintY + glintX + glintSize;
                    glintBackLeftIndex = width * (glintY + glintSize) + topLeftX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsRightEdge && !glintOverlapsTopEdge && !glintOverlapsBottomEdge) // right edge
                {
                    glintTopLeftIndex  = width * glintY + glintX;
                    glintTopRghtIndex  = width * glintY + backRightX;
                    glintBackLeftIndex = width * (glintY + glintSize) + glintX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                // check corner overlap
                
                if (glintOverlapsLeftEdge && glintOverlapsTopEdge) // top left corner
                {
                    glintTopLeftIndex  = topLeftIndex;
                    glintTopRghtIndex  = width * topLeftY + glintX + glintSize;
                    glintBackLeftIndex = width * (glintY + glintSize) + topLeftX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsRightEdge && glintOverlapsTopEdge) // top right corner
                {
                    glintTopLeftIndex  = width * topLeftY + glintX;
                    glintTopRghtIndex  = topRghtIndex ;
                    glintBackLeftIndex = width * (glintY + glintSize) + glintX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsLeftEdge && glintOverlapsBottomEdge) // bottom left corner
                {
                    glintTopLeftIndex  = width * glintY + topLeftX;
                    glintTopRghtIndex  = width * glintY + glintX + glintSize;
                    glintBackLeftIndex = backLeftIndex;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsRightEdge && glintOverlapsBottomEdge) // bottom right corner
                {
                    glintTopLeftIndex  = width * glintY + glintX;
                    glintTopRghtIndex  = width * glintY + backRightX;
                    glintBackLeftIndex = width * backRightY + glintX;
                    glintBackRghtIndex = backRghtIndex;
                }
                
                // calculate area and intensity of glint
                glintIntensity = I[glintBackRghtIndex] - I[glintBackLeftIndex] - I[glintTopRghtIndex] + I[glintTopLeftIndex];
                glintArea      = (glintTopRghtIndex - glintTopLeftIndex) * ((glintBackLeftIndex - glintTopLeftIndex) / width);
            }
            
            // calculate average pupil intensity, adjusting for glint
            
            double pupilIntensity = ((I[backRghtIndex] - I[backLeftIndex] - I[topRghtIndex] + I[topLeftIndex]) - glintIntensity) / (pupilArea - glintArea);
            
            if (pupilIntensity < minPupilIntensity)
            {
                pupilX = topLeftX;
                pupilY = topLeftY;
                minPupilIntensity = pupilIntensity;
            }
        }
    }
    
    haarProperties pupil;
    
    pupil.xPos = pupilX;
    pupil.yPos = pupilY;
    
    return pupil;
}

std::vector<int> radialGradient(const cv::Mat& img, int kernelSize, double pupilXCentre, double pupilYCentre)
{
    int kernelPerimeter = 8; // perimeter length of kernel
    double fc =  6.0;
    double sd =  1.0;

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

    uchar *ptr = img.data;
    int width  = img.cols;
    int height = img.rows;

    std::vector<int> gradientVector(width * height, 0);

    int borderSize = (kernelSize - 1) / 2;

    for (int y = borderSize; y < height - borderSize; y++)
    {
        for (int x = borderSize; x < width - borderSize; x++)
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

            double alpha = theta * kernelPerimeter / (2 * M_PI);
            double val = 0;

            for (int i = 0; i < kernelPerimeter; i++)
            {
                double dpos = std::abs(i - alpha);
                if (dpos > kernelPerimeter / 2) { dpos = kernelPerimeter - dpos; }
                double dneg = 0.5 * kernelPerimeter - dpos;

                double kernelVal = fc * (exp(- pow(dpos, 2) / sd) - exp(- pow(dneg, 2) / sd));

                // Do convolution
                int iPixel = width * (y + dY[i] * borderSize) + (x + dX[i] * borderSize);
                double pixelIntensity = ptr[iPixel];
                val += pixelIntensity * kernelVal;
            }

            if (val < 0) { val = 0; }
            int iPixel = width * y + x;
            gradientVector[iPixel] = val;
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

                        if (xOpposite < 0 || xOpposite >= haarWidth  || // ... and opposite pixel is out-of-bounds
                                yOpposite < 0 || yOpposite >= haarHeight ||
                                binaryImageVector[iOpposite] ==  0       || // ... or unfilled ...
                                binaryImageVector[iOpposite] == -1)
                        {
                            binaryImageVector[iCentre] = -1; // ... then remove pixel from edge
                        }
                    }
                }
            }
        }
    }

    return binaryImageVector;
}

std::vector<int> getEdgeIndices(const std::vector<int>& binaryImageVector)
{
    int haarSize = binaryImageVector.size();
    std::vector<int> cannyEdgeIndices;
    for (int i = 0; i < haarSize; i++)
    { if (binaryImageVector[i] == 1) { cannyEdgeIndices.push_back(i); }}
    return cannyEdgeIndices;
}

std::vector<edgeProperties> edgeFilter(const cv::Mat& img, std::vector<int>& p, int haarWidth, int haarHeight, double curvatureLowerLimit, double curvatureUpperLimit)
{
    std::vector<edgeProperties> vEdgePropertiesAll; // new structure containing length and indices of all edges
    
    int haarArea = haarWidth * haarHeight;
    
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

    int startEdgeSearchIndex = 0;
    
    while (true)
    {
        bool newEdgeFound = false;
        
        for (int i = startEdgeSearchIndex; i < haarArea; i++) // find a point belonging to an edge
        {
            if (p[i] == 1)
            {
                startEdgeSearchIndex = i;
                newEdgeFound = true;
                break;
            }
        }
        
        if (!newEdgeFound) { break; } // no (more) edges found
        
        // find edge endpoint
        
        std::vector<int> allEdgeIndicesRaw(1); // all indices belonging to the edge, out-of-order
        allEdgeIndicesRaw[0] = startEdgeSearchIndex;
        
        p[startEdgeSearchIndex] = 2; // tag pixel
        
        {
            std::vector<int> rawEdgeIndices(1);
            rawEdgeIndices[0] = startEdgeSearchIndex;
            
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
        
        int numberOfBranches = 1; // start off with one branch

        while (numberOfBranches > 0)
        {
            std::vector<int> branchSizes;                            // record length of each branch
            std::vector<std::vector<int>> branchIndicesVectors;      // record indices of each branch
            std::vector<std::vector<int>> branchStartIndicesVectors; // record start index of each branch

            for (int iBranch = 0; iBranch < numberOfBranches; iBranch++) // stops when junction is encountered
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
            numberOfBranches = branchStartIndices.size();
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
        
        // calculate curvature
        
        int edgeLength = allEdgeIndices.size();

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

            if (vectorAngle > M_PI)
            {
                vectorAngle = vectorAngle - 2 * M_PI;
            }
            else if (vectorAngle < -M_PI)
            {
                vectorAngle = vectorAngle + 2 * M_PI;
            }

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
                double edgePointCurvature = edgePointCurvatures[iEdgePoint];
                
                if (edgePointCurvature > 0)
                {
                    numberOfPositiveCurvatures++;
                }
                else if (edgePointCurvature < 0)
                {
                    numberOfNegativeCurvatures++;
                }
            }
            
            if (numberOfPositiveCurvatures >= numberOfNegativeCurvatures)
            {
                edgeCurvatureSign =  1;
            }
            else
            {
                edgeCurvatureSign = -1;
            }
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

            if (partialEdgeLength < curvatureWindowLength) { continue; } // ignore short edges

            // grab indices of edge points

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
                    {
                        partialEdgeIntensities[iEdgePoint] = ptr_img[edgePointXPos + edgePointYPos * haarWidth];
                    }
                    else
                    {
                        partialEdgeIntensities[iEdgePoint] = ptr_img[offsetXPos + offsetYPos * haarWidth];
                    }
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
            mEdgeProperties.length       = partialEdgeLength;  // add edge length to structure
            mEdgeProperties.size         = numEdgePoints;      // add total number of edge points of edge to structure
            mEdgeProperties.pointIndices = partialEdgeIndices;

            vEdgePropertiesAll.push_back(mEdgeProperties);
        }
    }

    return vEdgePropertiesAll; // return structure
}

std::vector<int> edgeThreshold(eyeProperties mEyeProperties, const std::vector<edgeProperties>& vEdgePropertiesAll)
{
    double expectedCurvature = mEyeProperties.v.edgeCurvaturePrediction;

    int numEdgesMax = mEyeProperties.p.edgeMaximumFitNumber;
    int numEdges    = vEdgePropertiesAll.size();

    std::vector<int> acceptedEdges(numEdgesMax);

    std::vector<double> totalScoresUnsorted(numEdges);

    for (int iEdge = 0; iEdge < numEdges; iEdge++)
    {
        // intensity score
        double intensity = vEdgePropertiesAll[iEdge].intensity;
        double intensityScore = (20 / (1 + 0.01 * pow(0.90, -intensity + mEyeProperties.v.edgeIntensityPrediction)));
        if (intensityScore < 0) { intensityScore = 0; }

        // position score
        double dR = std::abs(vEdgePropertiesAll[iEdge].distance - radiusPredictionFactor * mEyeProperties.v.radiusPrediction); // reduce prediction by constant factor, because actual pupil edge should envelop prediction
        double positionScore = (-15 / (mEyeProperties.v.radiusPrediction)) * dR + 15;
        if (positionScore < 0) { positionScore = 0; }

        // length score
        double length = vEdgePropertiesAll[iEdge].length;
        double lengthScore;
        if (length <= mEyeProperties.v.circumferencePrediction)
        {
            lengthScore = 12 * (1 - exp(-0.0002 * mEyeProperties.v.circumferencePrediction * length)); // longer = better
        }
        else
        {
            lengthScore = 12 / (1 + 0.01 * pow(0.85, -length + mEyeProperties.v.circumferencePrediction)); // closer to prediction = better
        }
        if (lengthScore < 0) { lengthScore = 0; }

        // curvature score
        double dC = std::abs(vEdgePropertiesAll[iEdge].curvatureAvg - expectedCurvature); // closer to prediction = better
        double curvatureScore = (-7 / expectedCurvature) * dC + 7;
        if (curvatureScore < 0) { curvatureScore = 0; }

        totalScoresUnsorted[iEdge] = intensityScore + positionScore + lengthScore + curvatureScore;
    }

    // Grab edges with highest score

    std::vector<double> totalScoresSorted = totalScoresUnsorted;
    std::sort(totalScoresSorted.begin(), totalScoresSorted.end());
    std::reverse(totalScoresSorted.begin(), totalScoresSorted.end());

    for (int iEdge = 0; iEdge < numEdgesMax; iEdge++) // THRESHOLD: do not handle more than a fixed number of edges
    {
        for (int jEdge = 0; jEdge < numEdges; jEdge++)
        {
            if (totalScoresSorted[iEdge] == totalScoresUnsorted[jEdge])
            {
                acceptedEdges[iEdge] = jEdge;

                break;
            }
        }
    }

    return acceptedEdges;

}

std::vector<double> EllipseRotationTransformation(const std::vector<double>& c)
{
    double A = c[0];
    double B = c[1];
    double C = c[2];
    double D = c[3];
    double E = c[4];
    double F = c[5];
    
    double alpha = 0.5 * (atan(B / (A - C))); // rotation angle
    
    double AA =  A * cos(alpha) * cos(alpha) + B * cos(alpha) * sin(alpha) + C * sin(alpha) * sin(alpha);
    double CC =  A * sin(alpha) * sin(alpha) - B * cos(alpha) * sin(alpha) + C * cos(alpha) * cos(alpha);
    double DD =  D * cos(alpha) + E * sin(alpha);
    double EE = -D * sin(alpha) + E * cos(alpha);
    double FF =  F;
    
    // semi axes
    
    double a = sqrt((-4 * FF * AA * CC + CC * DD * DD + AA * EE * EE)/(4 * AA * CC * CC));
    double b = sqrt((-4 * FF * AA * CC + CC * DD * DD + AA * EE * EE)/(4 * AA * AA * CC));
    
    // coordinates of centre point
    
    double x = -(DD / (2 * AA)) * cos(alpha) + (EE / (2 * CC)) * sin(alpha);
    double y = -(DD / (2 * AA)) * sin(alpha) - (EE / (2 * CC)) * cos(alpha);

    // width and height

    double w = 2 * sqrt(pow(a * cos(alpha), 2) + pow(b * sin(alpha), 2));
    double h = 2 * sqrt(pow(a * sin(alpha), 2) + pow(b * cos(alpha), 2));

    std::vector<double> v(6);
    
    v[0] = a;
    v[1] = b;
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

    ellipseProperties mEllipseProperties;
    mEllipseProperties.pupilDetected = true;

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

    Eigen::VectorXd EigenValues = EigenSolver.eigenvalues().real();
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

    double semiAxisA = ellipseParameters[0];
    double semiAxisB = ellipseParameters[1];

    double majorAxis = 0;
    double minorAxis = 0;

    if (semiAxisA >= semiAxisB)
    {
        majorAxis = semiAxisA;
        minorAxis = semiAxisB;
    }
    else
    {
        majorAxis = semiAxisB;
        minorAxis = semiAxisA;
    }

    if (!std::isnormal(majorAxis) || !std::isnormal(minorAxis)) // only take square root of positive values
    { mEllipseProperties.pupilDetected = false; }

    double h = pow((majorAxis - minorAxis), 2) / pow((majorAxis + minorAxis), 2);
    mEllipseProperties.circumference = M_PI * (majorAxis + minorAxis) * (1 + (3 * h) / (10 + sqrt(4 - 3 * h))); // ramanujans 2nd approximation
    mEllipseProperties.aspectRatio   = minorAxis / majorAxis;
    mEllipseProperties.radius        = 0.5 * (minorAxis + majorAxis);
    mEllipseProperties.xPos          = ellipseParameters[2];
    mEllipseProperties.yPos          = ellipseParameters[3];
    mEllipseProperties.width         = ellipseParameters[4];
    mEllipseProperties.height        = ellipseParameters[5];
    mEllipseProperties.coefficients  = ellipseFitCoefficients;

    return mEllipseProperties;
}

ellipseProperties findBestEllipseFit(const std::vector<edgeProperties>& vEdgePropertiesAll, int haarWidth, eyeProperties mEyeProperties, std::vector<int> cannyEdgePointIndices)
{
    ellipseProperties mEllipseProperties;
    mEllipseProperties.pupilDetected = false;
    
    int totalNumberOfEdges = vEdgePropertiesAll.size(); // total number of edges
    int totalNumberOfEdgePoints = cannyEdgePointIndices.size();

    std::vector<ellipseProperties> vEllipsePropertiesAll; // vector to record information for each accepted ellipse fit
    
    int numberOfFits = 0; // iterator over all accepted combinations

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
            std::vector<int>    combiEdgeIndices(combiNumEdges);
            std::vector<int>    combiEdgeLengths(combiNumEdges);
            std::vector<int>    combiEdgeSizes(combiNumEdges);
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
            
            // calculate lengths

            int edgeSetLength = std::accumulate(combiEdgeLengths.begin(), combiEdgeLengths.end(), 0);
            
            if (edgeSetLength < mEyeProperties.v.circumferencePrediction * edgeCollectionFraction)
            { continue; }// ignore edge collections that are too short

            int edgeSetSize = std::accumulate(combiEdgeSizes.begin(), combiEdgeSizes.end(), 0);

            // concatenate index vectors

            std::vector<int> edgeIndices; // vector containing all indices for fit
            edgeIndices.reserve(edgeSetSize); // preallocate memory

            for (int iEdge = 0; iEdge < combiNumEdges; iEdge++)
            { edgeIndices.insert(edgeIndices.end(), combiEdgePointIndices[iEdge].begin(), combiEdgePointIndices[iEdge].end()); }

            // fit ellipse

            ellipseProperties mEllipsePropertiesNew = fitEllipse(edgeIndices, edgeSetSize, haarWidth);

            if (!mEllipsePropertiesNew.pupilDetected) // error
            { continue; }

            // Size and shape filters
            
            if (mEllipsePropertiesNew.circumference > mEyeProperties.p.circumferenceMax || mEllipsePropertiesNew.circumference < mEyeProperties.p.circumferenceMin) // no large or small pupils
            { continue; }

            if (mEllipsePropertiesNew.aspectRatio < mEyeProperties.p.aspectRatioMin) // no extreme deviations from circular shape
            { continue; }

            if (std::abs(mEllipsePropertiesNew.circumference - mEyeProperties.v.circumferencePrediction) > mEyeProperties.v.thresholdCircumferenceChange) // no large pupil size changes
            { continue; }

            if (std::abs(mEllipsePropertiesNew.aspectRatio - mEyeProperties.v.aspectRatioPrediction) > mEyeProperties.v.thresholdAspectRatioChange) // no large pupil size changes
            { continue; }

            // calculate error between fit and every edge point

            std::vector<int> fittedEdgePoints;
            int fitSetLength;

            {
                double A = mEllipsePropertiesNew.coefficients[0];
                double B = mEllipsePropertiesNew.coefficients[1];
                double C = mEllipsePropertiesNew.coefficients[2];
                double D = mEllipsePropertiesNew.coefficients[3];
                double E = mEllipsePropertiesNew.coefficients[4];
                double F = mEllipsePropertiesNew.coefficients[5];

                // iterate over all edge points of all raw canny edges

                for (int iEdgePoint = 0; iEdgePoint < totalNumberOfEdgePoints; iEdgePoint++)
                {
                    int edgePointIndex = cannyEdgePointIndices[iEdgePoint];

                    double x = edgePointIndex % haarWidth;
                    double y = (edgePointIndex - x) / haarWidth;

                    double fitError = std::abs(A * x * x + B * x * y + C * y * y + D * x + E * y + F);

                    if (fitError < mEyeProperties.p.ellipseFitErrorMaximum)
                    { fittedEdgePoints.push_back(edgePointIndex); }
                }

                fitSetLength = fittedEdgePoints.size();

                if (fitSetLength < mEyeProperties.v.circumferencePrediction * minimumFitFraction)
                { continue; }// Threshold: ignore edge collections that are too short
            }

            // do refinement if ellipse fit overlaps with many new edge points

            if ((fitSetLength - edgeSetLength) > mEyeProperties.v.circumferencePrediction * minimumRefinementFraction)
            { mEllipsePropertiesNew = fitEllipse(fittedEdgePoints, fitSetLength, haarWidth); }

            // save parameters of accepted fit

            mEllipsePropertiesNew.edgeIndices = combiEdgeIndices;
            mEllipsePropertiesNew.intensity   = calculateMean(combiEdgeIntensities);
            mEllipsePropertiesNew.edgeLength  = fitSetLength;
            vEllipsePropertiesAll.push_back(mEllipsePropertiesNew);
            
            numberOfFits++;
        }
        while (std::next_permutation(edgeCombination.begin(), edgeCombination.end()));
    }
    
    // grab ellipse fit that resembles prior the most in size and shape
    
    if (numberOfFits > 0)
    {
        std::vector<double> featureChange(numberOfFits); // new type of fit error

        for (int iFit = 0; iFit < numberOfFits; iFit++)
        {
            double circumferenceChange = (std::abs(vEllipsePropertiesAll[iFit].circumference - mEyeProperties.v.circumferencePrediction));

            if (circumferenceChange < 0.5 * mEyeProperties.v.thresholdCircumferenceChange)
            {   circumferenceChange = 0.5 * mEyeProperties.v.thresholdCircumferenceChange; }

            double aspectRatioChange = (std::abs(vEllipsePropertiesAll[iFit].aspectRatio - mEyeProperties.v.aspectRatioPrediction));

            if (aspectRatioChange < 0.5 * mEyeProperties.v.thresholdAspectRatioChange)
            {   aspectRatioChange = 0.5 * mEyeProperties.v.thresholdAspectRatioChange; }

            double normalizationConstant = mEyeProperties.v.thresholdCircumferenceChange / mEyeProperties.v.thresholdAspectRatioChange;
            aspectRatioChange = normalizationConstant * aspectRatioChange;

            featureChange[iFit] = (circumferenceChange + aspectRatioChange) / sqrt((double) vEllipsePropertiesAll[iFit].edgeLength);
        }

        int acceptedFitIndex = std::distance(featureChange.begin(), std::min_element(featureChange.begin(), featureChange.end()));

        mEllipseProperties = vEllipsePropertiesAll[acceptedFitIndex];
        mEllipseProperties.pupilDetected = true;
    }
    
    return mEllipseProperties;
}

eyeProperties pupilDetection(const cv::Mat& imageOriginalBGR, eyeProperties mEyeProperties)
{
    // Define some variables
    
    eyeProperties mEyePropertiesNew = mEyeProperties; // new properties for new frame
    std::vector<edgeProperties> vEdgePropertiesNew;
    ellipseProperties mEllipseProperties;

    mEllipseProperties.pupilDetected  = false;
    mEyePropertiesNew.m.errorDetected = false;
    mEyePropertiesNew.m.image = imageOriginalBGR;

    // Define search area

    int imageWdth = imageOriginalBGR.cols;
    int imageHght = imageOriginalBGR.rows;

    int searchStartX = round(mEyeProperties.v.xPosPredicted - mEyeProperties.v.searchRadius);
    int searchStartY = round(mEyeProperties.v.yPosPredicted - mEyeProperties.v.searchRadius);
    
    if (searchStartX < 0) { searchStartX = 0; }
    if (searchStartY < 0) { searchStartY = 0; }
    
    int searchEndX = round(mEyeProperties.v.xPosPredicted + mEyeProperties.v.searchRadius);
    int searchEndY = round(mEyeProperties.v.yPosPredicted + mEyeProperties.v.searchRadius);
    
    if (searchEndX >= imageWdth) { searchEndX = imageWdth - 1; }
    if (searchEndY >= imageHght) { searchEndY = imageHght - 1; }
    
    int searchWdth = searchEndX - searchStartX + 1;
    int searchHght = searchEndY - searchStartY + 1;
    
    int pupilHaarWdth = round(pupilHaarReductionFactor * mEyeProperties.v.widthPrediction);
    int pupilHaarHght = round(pupilHaarReductionFactor * mEyeProperties.v.heightPrediction);

    int offsetPupilHaarXPos = 0;
    int offsetPupilHaarYPos = 0;
    int offsetPupilHaarWdth = 0;
    int offsetPupilHaarHght = 0;

    if (pupilHaarWdth > searchWdth) { pupilHaarWdth = searchWdth; }
    if (pupilHaarHght > searchHght) { pupilHaarHght = searchHght; }

    if (pupilHaarWdth > 0 && pupilHaarHght > 0)
    {
        // Needed for offline mode when threshold is too low

        if (mEyeProperties.v.thresholdCircumferenceChange > mEyeProperties.p.circumferenceMax)
        {   mEyeProperties.v.thresholdCircumferenceChange = mEyeProperties.p.circumferenceMax; }
        else if (mEyeProperties.v.thresholdCircumferenceChange < mEyeProperties.p.circumferenceChangeThreshold)
        {        mEyeProperties.v.thresholdCircumferenceChange = mEyeProperties.p.circumferenceChangeThreshold; }

        if (mEyeProperties.v.thresholdAspectRatioChange > 1.0)
        {   mEyeProperties.v.thresholdAspectRatioChange = 1.0; }
        else if (mEyeProperties.v.thresholdAspectRatioChange < mEyeProperties.p.aspectRatioChangeThreshold)
        {        mEyeProperties.v.thresholdAspectRatioChange = mEyeProperties.p.aspectRatioChangeThreshold; }

        // Convert to grayscale

        cv::Mat imageOriginalGray;
        cv::cvtColor(imageOriginalBGR, imageOriginalGray, cv::COLOR_BGR2GRAY);

        // Initial approximate detection

        std::vector<unsigned int> integralImage = calculateIntImg(imageOriginalGray, imageWdth, searchStartX, searchStartY, searchWdth, searchHght);

        int glintSize = mEyeProperties.p.glintSize;

        haarProperties glintHaarProperties = detectGlint(imageOriginalGray, imageWdth, searchStartX, searchStartY, searchWdth, searchHght, glintSize);

        int glintXPos = searchStartX + glintHaarProperties.xPos;
        int glintYPos = searchStartY + glintHaarProperties.yPos;

        haarProperties pupilHaarProperties = detectPupilApprox(integralImage, searchWdth, searchHght, pupilHaarWdth, pupilHaarHght, glintXPos, glintYPos, mEyeProperties.p.glintSize);

        int pupilHaarXPos = searchStartX + pupilHaarProperties.xPos;
        int pupilHaarYPos = searchStartY + pupilHaarProperties.yPos;

        offsetPupilHaarXPos = pupilHaarXPos -     mEyeProperties.p.pupilOffset;
        offsetPupilHaarYPos = pupilHaarYPos -     mEyeProperties.p.pupilOffset;
        offsetPupilHaarWdth = pupilHaarWdth + 2 * mEyeProperties.p.pupilOffset;
        offsetPupilHaarHght = pupilHaarHght + 2 * mEyeProperties.p.pupilOffset;

        // Check limits

        if (offsetPupilHaarWdth >= imageWdth) { offsetPupilHaarWdth = imageWdth - 1; }
        if (offsetPupilHaarHght >= imageHght) { offsetPupilHaarHght = imageHght - 1; }
        if (offsetPupilHaarXPos < 0) { offsetPupilHaarXPos = 0; }
        if (offsetPupilHaarYPos < 0) { offsetPupilHaarYPos = 0; }
        if (offsetPupilHaarXPos + offsetPupilHaarWdth >= imageWdth) { offsetPupilHaarWdth = imageWdth - offsetPupilHaarXPos - 1; }
        if (offsetPupilHaarYPos + offsetPupilHaarHght >= imageHght) { offsetPupilHaarHght = imageHght - offsetPupilHaarYPos - 1; }

        // Crop image to outer pupil Haar region

        cv::Rect pupilROI(offsetPupilHaarXPos, offsetPupilHaarYPos, offsetPupilHaarWdth, offsetPupilHaarHght);
        cv::Mat imagePupilBGR = imageOriginalBGR(pupilROI);

        // Convert back to grayscale

        cv::Mat imagePupilGray;
        cv::cvtColor(imagePupilBGR, imagePupilGray, cv::COLOR_BGR2GRAY);

        // Canny edge detection

        cv::Mat imagePupilGrayBlurred;
        int cannyBlurLevel = 2 * mEyeProperties.p.cannyBlurLevel - 1; // should be odd
        if (cannyBlurLevel > 0) { cv::GaussianBlur(imagePupilGray, imagePupilGrayBlurred, cv::Size(cannyBlurLevel, cannyBlurLevel), 0, 0);
        } else                  { imagePupilGrayBlurred = imagePupilGray; }

        double xPosPredictedRelative = mEyeProperties.v.xPosPredicted - offsetPupilHaarXPos;
        double yPosPredictedRelative = mEyeProperties.v.yPosPredicted - offsetPupilHaarYPos;

        std::vector<int>  imgGradient           = radialGradient(imagePupilGrayBlurred, mEyeProperties.p.cannyKernelSize, xPosPredictedRelative, yPosPredictedRelative);
        std::vector<int>  imgGradientSuppressed = nonMaximumSuppresion(imgGradient, offsetPupilHaarWdth, offsetPupilHaarHght, xPosPredictedRelative, yPosPredictedRelative);
        std::vector<int>  cannyEdges            = hysteresisTracking(imgGradientSuppressed, offsetPupilHaarWdth, offsetPupilHaarHght, mEyeProperties.p.cannyThresholdHigh, mEyeProperties.p.cannyThresholdLow);
        std::vector<int>  cannyEdgesSharpened   = sharpenEdges(cannyEdges, offsetPupilHaarWdth, offsetPupilHaarHght);
        std::vector<int>  cannyEdgesIndices     = getEdgeIndices(cannyEdgesSharpened);

        // Edge thresholding

        double curvatureLowerLimit;

        {
            double A =  5796;
            double B = -0.7855;
            double C = -0.7259;
            double D =  0.3801;
            double E = -1.408;

            curvatureLowerLimit = A * exp(B * pow(mEyeProperties.v.circumferencePrediction, D) + C * pow(mEyeProperties.v.aspectRatioPrediction, E)) - mEyeProperties.v.curvatureOffset;
        }

        double curvatureUpperLimit;

        {
            double A =  98.31;
            double B =  9.95;
            double C = -2.106;
            double D = -0.412;
            double E =  0.4122;

            curvatureUpperLimit = A * exp(B * pow(mEyeProperties.v.circumferencePrediction, D) + C * pow(mEyeProperties.v.aspectRatioPrediction, E)) + mEyeProperties.v.curvatureOffset;
        }

        std::vector<edgeProperties> vEdgePropertiesAll = edgeFilter(imagePupilGray, cannyEdgesSharpened, offsetPupilHaarWdth, offsetPupilHaarHght, curvatureLowerLimit, curvatureUpperLimit);

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
                    int edgePointXPos  =  edgePointIndex % offsetPupilHaarWdth;
                    int edgePointYPos  = (edgePointIndex - edgePointXPos) / offsetPupilHaarWdth;

                    double dX = mEyeProperties.v.xPosPredicted - (edgePointXPos + offsetPupilHaarXPos);
                    double dY = mEyeProperties.v.yPosPredicted - (edgePointYPos + offsetPupilHaarYPos);

                    dR += sqrt(pow(dX, 2) + pow(dY, 2));
                }

                vEdgePropertiesAll[iEdge].distance = dR / edgeSize;
            }

            mEyeProperties.v.edgeCurvaturePrediction = 0.5 * (curvatureUpperLimit + curvatureLowerLimit);

            std::vector<int> acceptedEdges = edgeThreshold(mEyeProperties, vEdgePropertiesAll);

            for (int iEdge = 0; iEdge < mEyeProperties.p.edgeMaximumFitNumber; iEdge++) // grab accepted edges
            {
                int jEdge = acceptedEdges[iEdge];
                vEdgePropertiesAll[jEdge].flag = 1;
                vEdgePropertiesNew.push_back(vEdgePropertiesAll[jEdge]);
            }

            mEllipseProperties = findBestEllipseFit(vEdgePropertiesNew, offsetPupilHaarWdth, mEyeProperties, cannyEdgesIndices); // ellipse fitting
        }

        // Classify edges

        int numEdgesNew = mEllipseProperties.edgeIndices.size();

        for (int iEdge = 0; iEdge < numEdgesNew; iEdge++)
        {
            int jEdge = mEllipseProperties.edgeIndices[iEdge];
            vEdgePropertiesAll[jEdge].flag = 2;
        }

        // Features for all edges

        mEyePropertiesNew.m.edgePropertiesAll = vEdgePropertiesAll;

        // Save parameters

        mEyePropertiesNew.v.edgeCurvaturePrediction = mEyeProperties.v.edgeCurvaturePrediction;

        mEyePropertiesNew.v.aspectRatioExact   = mEllipseProperties.aspectRatio;
        mEyePropertiesNew.v.circumferenceExact = mEllipseProperties.circumference;
        mEyePropertiesNew.v.pupilDetected      = mEllipseProperties.pupilDetected;

        // For draw functions

        mEyePropertiesNew.m.offsetPupilHaarXPos = offsetPupilHaarXPos;
        mEyePropertiesNew.m.offsetPupilHaarYPos = offsetPupilHaarYPos;
        mEyePropertiesNew.m.offsetPupilHaarWdth = offsetPupilHaarWdth;
        mEyePropertiesNew.m.offsetPupilHaarHght = offsetPupilHaarHght;

        mEyePropertiesNew.m.pupilHaarXPos = pupilHaarXPos;
        mEyePropertiesNew.m.pupilHaarYPos = pupilHaarYPos;
        mEyePropertiesNew.m.pupilHaarWdth = pupilHaarWdth;

        mEyePropertiesNew.m.glintXPos = glintXPos;
        mEyePropertiesNew.m.glintYPos = glintYPos;
        mEyePropertiesNew.m.glintSize = glintSize;

        mEyePropertiesNew.m.cannyEdges          = cannyEdges;
        mEyePropertiesNew.m.ellipseCoefficients = mEllipseProperties.coefficients;
    }
    else
    {
        mEyePropertiesNew.m.errorDetected = true;
    }

    // For running averages
    
    if (!mEllipseProperties.pupilDetected) // pupil not detected
    {
        mEyePropertiesNew.v.aspectRatioAverage    = mEyeProperties.v.aspectRatioAverage + mEyeProperties.p.alphaAverage * (mEyeProperties.v.aspectRatioPrediction - mEyeProperties.v.aspectRatioAverage);
        mEyePropertiesNew.v.aspectRatioMomentum   = mEyeProperties.v.aspectRatioMomentum * mEyeProperties.p.alphaMomentum;
        mEyePropertiesNew.v.aspectRatioPrediction = mEyeProperties.v.aspectRatioPrediction + mEyeProperties.p.alphaPrediction * (mEyePropertiesNew.v.aspectRatioAverage - mEyeProperties.v.aspectRatioPrediction);

        mEyePropertiesNew.v.circumferenceAverage    = mEyeProperties.v.circumferenceAverage + mEyeProperties.p.alphaAverage * (mEyeProperties.v.circumferencePrediction - mEyeProperties.v.circumferenceAverage);
        mEyePropertiesNew.v.circumferenceMomentum   = mEyeProperties.v.circumferenceMomentum * mEyeProperties.p.alphaMomentum;
        mEyePropertiesNew.v.circumferencePrediction = mEyeProperties.v.circumferencePrediction + mEyeProperties.p.alphaPrediction * (mEyePropertiesNew.v.circumferenceAverage - mEyeProperties.v.circumferencePrediction);

        mEyePropertiesNew.v.curvatureOffset = mEyeProperties.v.curvatureOffset * (2 - mEyeProperties.p.alphaMiscellaneous);

        mEyePropertiesNew.v.edgeIntensityAverage    = mEyeProperties.v.edgeIntensityAverage    + mEyeProperties.p.alphaAverage    * (mEyeProperties.v.edgeIntensityPrediction - mEyeProperties.v.edgeIntensityAverage);
        mEyePropertiesNew.v.edgeIntensityPrediction = mEyeProperties.v.edgeIntensityPrediction + mEyeProperties.p.alphaPrediction * (mEyeProperties.v.edgeIntensityAverage - mEyeProperties.v.edgeIntensityPrediction);

        mEyePropertiesNew.v.heightAverage    = mEyeProperties.v.heightAverage + mEyeProperties.p.alphaAverage * (mEyeProperties.v.heightPrediction - mEyeProperties.v.heightAverage);
        mEyePropertiesNew.v.heightMomentum   = mEyeProperties.v.heightMomentum * mEyeProperties.p.alphaMomentum;
        mEyePropertiesNew.v.heightPrediction = mEyeProperties.v.heightPrediction + mEyeProperties.p.alphaPrediction * (mEyePropertiesNew.v.heightAverage - mEyeProperties.v.heightPrediction);

        mEyePropertiesNew.v.radiusMomentum   = mEyeProperties.v.radiusMomentum * mEyeProperties.p.alphaMomentum;
        mEyePropertiesNew.v.radiusPrediction = mEyeProperties.v.circumferencePrediction / (2 * M_PI);

        mEyePropertiesNew.v.searchRadius = mEyeProperties.v.searchRadius * (2 - mEyeProperties.p.alphaMiscellaneous);

        mEyePropertiesNew.v.thresholdAspectRatioChange   = mEyeProperties.v.thresholdAspectRatioChange   * (2 - mEyeProperties.p.alphaMiscellaneous);
        mEyePropertiesNew.v.thresholdCircumferenceChange = mEyeProperties.v.thresholdCircumferenceChange * (2 - mEyeProperties.p.alphaMiscellaneous);

        mEyePropertiesNew.v.widthAverage    = mEyeProperties.v.widthAverage + mEyeProperties.p.alphaAverage * (mEyeProperties.v.widthPrediction - mEyeProperties.v.widthAverage);
        mEyePropertiesNew.v.widthMomentum   = mEyeProperties.v.widthMomentum * mEyeProperties.p.alphaMomentum;
        mEyePropertiesNew.v.widthPrediction = mEyeProperties.v.widthPrediction + mEyeProperties.p.alphaPrediction * (mEyePropertiesNew.v.widthAverage - mEyeProperties.v.widthPrediction);

        mEyePropertiesNew.v.xPosPredicted = mEyeProperties.v.xPosPredicted + mEyeProperties.p.alphaPrediction * (offsetPupilHaarXPos + 0.5 * offsetPupilHaarWdth - mEyeProperties.v.xPosPredicted) + mEyeProperties.v.xVelocity;
        mEyePropertiesNew.v.xVelocity     = mEyeProperties.v.xVelocity * mEyeProperties.p.alphaMomentum;

        mEyePropertiesNew.v.yPosPredicted = mEyeProperties.v.yPosPredicted + mEyeProperties.p.alphaPrediction * (offsetPupilHaarYPos + 0.5 * offsetPupilHaarHght - mEyeProperties.v.yPosPredicted) + mEyeProperties.v.yVelocity;
        mEyePropertiesNew.v.yVelocity     = mEyeProperties.v.yVelocity * mEyeProperties.p.alphaMomentum;
    }
    else // pupil detected
    {
        mEyePropertiesNew.v.aspectRatioAverage    =  mEyeProperties.v.aspectRatioAverage    + mEyeProperties.p.alphaAverage    * (mEyeProperties.v.aspectRatioPrediction - mEyeProperties.v.aspectRatioAverage);
        mEyePropertiesNew.v.aspectRatioMomentum   = (mEyeProperties.v.aspectRatioMomentum   + mEyeProperties.p.alphaMomentum   * (mEyePropertiesNew.v.aspectRatioPrediction - mEyeProperties.v.aspectRatioPrediction));
        mEyePropertiesNew.v.aspectRatioPrediction =  mEyeProperties.v.aspectRatioPrediction + mEyeProperties.p.alphaPrediction * (mEllipseProperties.aspectRatio - mEyeProperties.v.aspectRatioPrediction) + mEyeProperties.v.aspectRatioMomentum;

        mEyePropertiesNew.v.circumferenceAverage    =  mEyeProperties.v.circumferenceAverage    + mEyeProperties.p.alphaAverage    * (mEyeProperties.v.circumferencePrediction - mEyeProperties.v.circumferenceAverage);
        mEyePropertiesNew.v.circumferenceMomentum   = (mEyeProperties.v.circumferenceMomentum   + mEyeProperties.p.alphaMomentum   * (mEyePropertiesNew.v.circumferencePrediction - mEyeProperties.v.circumferencePrediction));
        mEyePropertiesNew.v.circumferencePrediction =  mEyeProperties.v.circumferencePrediction + mEyeProperties.p.alphaPrediction * (mEllipseProperties.circumference - mEyeProperties.v.circumferencePrediction) + mEyeProperties.v.circumferenceMomentum;

        mEyePropertiesNew.v.curvatureOffset = mEyeProperties.v.curvatureOffset * mEyeProperties.p.alphaMiscellaneous;

        mEyePropertiesNew.v.edgeIntensityAverage    = mEyeProperties.v.edgeIntensityAverage    + mEyeProperties.p.alphaAverage    * (mEyeProperties.v.edgeIntensityPrediction - mEyeProperties.v.edgeIntensityAverage);
        mEyePropertiesNew.v.edgeIntensityPrediction = mEyeProperties.v.edgeIntensityPrediction + mEyeProperties.p.alphaPrediction * (mEllipseProperties.intensity - mEyeProperties.v.edgeIntensityPrediction);

        mEyePropertiesNew.v.heightAverage    =  mEyeProperties.v.heightAverage    + mEyeProperties.p.alphaAverage    * (mEyeProperties.v.heightPrediction - mEyeProperties.v.heightAverage);
        mEyePropertiesNew.v.heightMomentum   = (mEyeProperties.v.heightMomentum   + mEyeProperties.p.alphaMomentum   * (mEyePropertiesNew.v.heightPrediction - mEyeProperties.v.heightPrediction));
        mEyePropertiesNew.v.heightPrediction =  mEyeProperties.v.heightPrediction + mEyeProperties.p.alphaPrediction * (mEllipseProperties.height - mEyeProperties.v.heightPrediction) + mEyeProperties.v.heightMomentum;

        mEyePropertiesNew.v.radiusMomentum   = (mEyeProperties.v.radiusMomentum + (mEyePropertiesNew.v.radiusPrediction - mEyeProperties.v.radiusPrediction)) * mEyeProperties.p.alphaMomentum;
        mEyePropertiesNew.v.radiusPrediction =  mEyeProperties.v.radiusPrediction + mEyeProperties.p.alphaPrediction * (mEllipseProperties.radius - mEyeProperties.v.radiusPrediction) + mEyeProperties.v.radiusMomentum;

        mEyePropertiesNew.v.searchRadius = mEyeProperties.v.searchRadius * mEyeProperties.p.alphaMiscellaneous;

        mEyePropertiesNew.v.thresholdAspectRatioChange   = mEyeProperties.v.thresholdAspectRatioChange   * mEyeProperties.p.alphaMiscellaneous;
        mEyePropertiesNew.v.thresholdCircumferenceChange = mEyeProperties.v.thresholdCircumferenceChange * mEyeProperties.p.alphaMiscellaneous;

        mEyePropertiesNew.v.widthAverage    =  mEyeProperties.v.widthAverage    + mEyeProperties.p.alphaAverage    * (mEyeProperties.v.widthPrediction - mEyeProperties.v.widthAverage);
        mEyePropertiesNew.v.widthMomentum   = (mEyeProperties.v.widthMomentum   + mEyeProperties.p.alphaMomentum   * (mEyePropertiesNew.v.widthPrediction - mEyeProperties.v.widthPrediction));
        mEyePropertiesNew.v.widthPrediction =  mEyeProperties.v.widthPrediction + mEyeProperties.p.alphaPrediction * (mEllipseProperties.width - mEyeProperties.v.widthPrediction) + mEyeProperties.v.widthMomentum;

        mEyePropertiesNew.v.xPosExact = mEllipseProperties.xPos + offsetPupilHaarXPos;
        mEyePropertiesNew.v.xPosPredicted = mEyeProperties.v.xPosPredicted + mEyeProperties.p.alphaPrediction * (mEyePropertiesNew.v.xPosExact - mEyeProperties.v.xPosPredicted) + mEyeProperties.v.xVelocity;
        mEyePropertiesNew.v.xVelocity = (mEyeProperties.v.xVelocity + (mEyePropertiesNew.v.xPosPredicted - mEyeProperties.v.xPosPredicted)) * mEyeProperties.p.alphaMomentum;

        mEyePropertiesNew.v.yPosExact = mEllipseProperties.yPos + offsetPupilHaarYPos;
        mEyePropertiesNew.v.yPosPredicted = mEyeProperties.v.yPosPredicted + mEyeProperties.p.alphaPrediction * (mEyePropertiesNew.v.yPosExact - mEyeProperties.v.yPosPredicted) + mEyeProperties.v.yVelocity;
        mEyePropertiesNew.v.yVelocity = (mEyeProperties.v.yVelocity + (mEyePropertiesNew.v.yPosPredicted - mEyeProperties.v.yPosPredicted)) * mEyeProperties.p.alphaMomentum;

        // Grab pupil image

        int radius = round(pupilImageFactor * mEllipseProperties.radius);
        int xPos = mEyePropertiesNew.v.xPosExact - radius;
        int yPos = mEyePropertiesNew.v.yPosExact - radius;
        int wdth = 2 * radius;
        int hght = wdth;

        if (xPos < 0) { xPos = 0; }
        if (yPos < 0) { yPos = 0; }
        if (xPos + wdth >= imageWdth) { wdth = imageWdth - xPos - 1; }
        if (yPos + hght >= imageHght) { hght = imageHght - yPos - 1; }

        cv::Rect pupilROI(xPos, yPos, wdth, hght);
        cv::Mat imagePupilBGR = imageOriginalBGR(pupilROI);
        cv::Mat imagePupilGray;
        cv::cvtColor(imagePupilBGR, imagePupilGray, cv::COLOR_BGR2GRAY);

        mEyePropertiesNew.m.imagePupil = imagePupilGray;
    }
    
    // Check variable limits

    if (mEyePropertiesNew.v.searchRadius > imageWdth)
    {   mEyePropertiesNew.v.searchRadius = imageWdth; }
    else if (mEyePropertiesNew.v.searchRadius < (0.5 * offsetPupilHaarWdth)) // add height as well
    {        mEyePropertiesNew.v.searchRadius = ceil(0.5 * offsetPupilHaarWdth); }

    if (mEyePropertiesNew.v.thresholdCircumferenceChange > mEyePropertiesNew.p.circumferenceMax)
    {   mEyePropertiesNew.v.thresholdCircumferenceChange = mEyePropertiesNew.p.circumferenceMax; }
    else if (mEyePropertiesNew.v.thresholdCircumferenceChange < mEyePropertiesNew.p.circumferenceChangeThreshold)
    {        mEyePropertiesNew.v.thresholdCircumferenceChange = mEyePropertiesNew.p.circumferenceChangeThreshold; }

    if (mEyePropertiesNew.v.thresholdAspectRatioChange > 1.0) {
        mEyePropertiesNew.v.thresholdAspectRatioChange = 1.0; }
    else if (mEyePropertiesNew.v.thresholdAspectRatioChange < mEyePropertiesNew.p.aspectRatioChangeThreshold)
    {        mEyePropertiesNew.v.thresholdAspectRatioChange = mEyePropertiesNew.p.aspectRatioChangeThreshold; }

    if (mEyePropertiesNew.v.curvatureOffset > 180)
    {   mEyePropertiesNew.v.curvatureOffset = 180; }
    else if (mEyePropertiesNew.v.curvatureOffset < mEyePropertiesNew.p.curvatureOffsetMin)
    {        mEyePropertiesNew.v.curvatureOffset = mEyePropertiesNew.p.curvatureOffsetMin; }

    return mEyePropertiesNew;
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
