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
    if (value < 0.0)
    {
        return floor(value);
    }
    else
    {
        return ceil(value);
    }
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
            
            if (x == 0 && y == 0) // first point
            {
                integralImage[i] = val;
            }
            else if (y == 0) // first row
            {
                integralImage[i] = val + integralImage[i - 1];
            }
            else if (x == 0) // first column
            {
                integralImage[i] = val + integralImage[i - width];
            }
            else
            {
                integralImage[i] = val + integralImage[i - 1] + integralImage[i - width] - integralImage[i - width - 1];
            }
        }
    }
    
    return integralImage;
}

haarProperties GlintDetector(const cv::Mat& img, int imgWidth, int startX, int startY, int width, int height, int gradientWindowLength, int glintRadius)
{
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
    
    glint.x_pos = x - glintRadius;
    glint.y_pos = y - glintRadius;
    
    return glint;
}

haarProperties PupilHaarDetector(const std::vector<unsigned int>& I, int width, int height, int haarWidth, int glintX, int glintY, int glintRadius)
{
    int glintDiameter = 2 * glintRadius;
    
    int pupilX = 0;
    int pupilY = 0;
    
    double minPupilIntensity = pow(10, 10); // arbitrarily large value
    
    int pupilArea = (haarWidth - 1) * (haarWidth - 1);
    
    for (int y = 0; y < height - haarWidth; y++)
    {
        for (int x = 0; x < width - haarWidth; x++)
        {
            // vertices of inner square
            
            int topLeftX = x;
            int topLeftY = y;
            
            int backRightX = topLeftX + (haarWidth - 1);
            int backRightY = topLeftY + (haarWidth - 1);
            
            int topLeftIndex = width * topLeftY + topLeftX;
            int topRghtIndex = topLeftIndex + (haarWidth - 1);
            int backLeftIndex = topLeftIndex + (haarWidth - 1) * width;
            int backRghtIndex = topRghtIndex + backLeftIndex - topLeftIndex;
            
            // calculate glint intensity
            
            double glintIntensity = 0.0;
            double glintArea = 0.0;
            
            bool glintWithinHaarDetector = false; // flag for glint overlap
            
            std::vector<int> z(4);
            z[0] = 0;
            z[1] = 0;
            z[2] = glintDiameter;
            z[3] = glintDiameter;
            
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
                bool glintOverlapsLeftEdge = false;
                bool glintOverlapsRightEdge = false;
                bool glintOverlapsTopEdge = false;
                bool glintOverlapsBottomEdge = false;
                
                if (glintX < topLeftX)
                {
                    glintOverlapsLeftEdge = true;
                }
                else if (glintX + glintDiameter > backRightX)
                {
                    glintOverlapsRightEdge = true;
                }
                
                if (glintY < topLeftY)
                {
                    glintOverlapsTopEdge = true;
                }
                else if (glintY + glintDiameter > backRightY)
                {
                    glintOverlapsBottomEdge = true;
                }
                
                /// coordinates of corners of glint square
                int glintTopLeftIndex = width * glintY + glintX;
                int glintTopRghtIndex = width * glintY + glintX + glintDiameter;
                int glintBackLeftIndex = width * (glintY + glintDiameter) + glintX;
                int glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                
                // check if glint square overlaps with edge or corner of pupil square
                
                // check edge overlap
                
                if (!glintOverlapsLeftEdge && !glintOverlapsRightEdge && glintOverlapsTopEdge) // top edge
                {
                    glintTopLeftIndex = width * topLeftY + glintX;
                    glintTopRghtIndex = width * topLeftY + glintX + glintDiameter;
                    glintBackLeftIndex = width * (glintY + glintDiameter) + glintX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                
                if (!glintOverlapsLeftEdge && !glintOverlapsRightEdge && glintOverlapsBottomEdge) // bottom edge
                {
                    glintTopLeftIndex = width * glintY + glintX;
                    glintTopRghtIndex = width * glintY + glintX + glintDiameter;
                    glintBackLeftIndex = width * backRightY + glintX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                
                if (glintOverlapsLeftEdge && !glintOverlapsTopEdge && !glintOverlapsBottomEdge) // left edge
                {
                    glintTopLeftIndex = width * glintY + topLeftX;
                    glintTopRghtIndex = width * glintY + glintX + glintDiameter;
                    glintBackLeftIndex = width * (glintY + glintDiameter) + topLeftX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsRightEdge && !glintOverlapsTopEdge && !glintOverlapsBottomEdge) // right edge
                {
                    glintTopLeftIndex = width * glintY + glintX;
                    glintTopRghtIndex = width * glintY + backRightX;
                    glintBackLeftIndex = width * (glintY + glintDiameter) + glintX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                // check corner overlap
                
                if (glintOverlapsLeftEdge && glintOverlapsTopEdge) // top left corner
                {
                    glintTopLeftIndex = topLeftIndex;
                    glintTopRghtIndex = width * topLeftY + glintX + glintDiameter;
                    glintBackLeftIndex = width * (glintY + glintDiameter) + topLeftX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsRightEdge && glintOverlapsTopEdge) // top right corner
                {
                    glintTopLeftIndex = width * topLeftY + glintX;
                    glintTopRghtIndex = topRghtIndex ;
                    glintBackLeftIndex = width * (glintY + glintDiameter) + glintX;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsLeftEdge && glintOverlapsBottomEdge) // bottom left corner
                {
                    glintTopLeftIndex = width * glintY + topLeftX;
                    glintTopRghtIndex = width * glintY + glintX + glintDiameter;
                    glintBackLeftIndex = backLeftIndex;
                    glintBackRghtIndex = glintTopRghtIndex + glintBackLeftIndex - glintTopLeftIndex;
                }
                
                if (glintOverlapsRightEdge && glintOverlapsBottomEdge) // bottom right corner
                {
                    glintTopLeftIndex = width * glintY + glintX;
                    glintTopRghtIndex = width * glintY + backRightX;
                    glintBackLeftIndex = width * backRightY + glintX;
                    glintBackRghtIndex = backRghtIndex;
                }
                
                // calculate area and intensity of glint
                glintIntensity = I[glintBackRghtIndex] - I[glintBackLeftIndex] - I[glintTopRghtIndex] + I[glintTopLeftIndex];
                glintArea = (glintTopRghtIndex - glintTopLeftIndex) * ((glintBackLeftIndex - glintTopLeftIndex) / width);
            }
            
            /// calculate average pupil intensity, adjusting for glint
            
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
    
    pupil.x_pos = pupilX;
    pupil.y_pos = pupilY;
    
    return pupil;
}

std::vector<char> cannyConversion(const cv::Mat& img, int haarWidth)
{
    int haarSize = haarWidth * haarWidth;
    
    uchar *ptr_img = img.data;
    
    std::vector<char> binaryImageVectorRaw(haarSize);
    
    for (int i = 0; i < haarSize; i++)
    {
        if (ptr_img[i] == 255)
        {
            binaryImageVectorRaw[i] = 1;
        }
        else
        {
            binaryImageVectorRaw[i] = 0;
        }
    }
    
    return binaryImageVectorRaw;
}

std::vector<char> MorphOpen(std::vector<char>& binaryImageVectorRaw, int haarWidth)
{
    std::vector<char> binaryImageVector = binaryImageVectorRaw;

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

    for (int yCentre = 0; yCentre < haarWidth; yCentre++)
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

                    for (char n = 0; n < 2; n++) // loop through two connected neighbouring pixels
                    {
                        char q = 2 * (m + n) % 8;

                        int xNeighbour = xCentre + dX[q];
                        int yNeighbour = yCentre + dY[q];

                        if (xNeighbour < 0 || xNeighbour >= haarWidth ||yNeighbour < 0 || yNeighbour >= haarWidth)
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

                        if (xOpposite < 0 || xOpposite >= haarWidth || // ... and opposite pixel is out-of-bounds
                                yOpposite < 0 || yOpposite >= haarWidth ||
                                binaryImageVector[iOpposite] ==  0 || // ... or unfilled ...
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

std::vector<int> findCannyIndices(const std::vector<char>& binaryImageVector)
{
    int haarSize = binaryImageVector.size();

    std::vector<int> cannyEdgeIndices;

    for (int i = 0; i < haarSize; i++)
    {
        if (binaryImageVector[i] == 1)
        {
            cannyEdgeIndices.push_back(i);
        }
    }

    return cannyEdgeIndices;
}

edgeProperties EdgeFilter(const cv::Mat& img, std::vector<char>& p, int haarWidth, double curvatureLowerLimit, double curvatureUpperLimit)
{
    edgeProperties mEdgePropertiesAll; // new structure containing length and indices of all edges
    
    int haarArea = haarWidth * haarWidth;
    
    /// scanning neighbours
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
        
        if (!newEdgeFound)
        {
            break; // no (more) edges found
        }
        
        // find edge end point
        
        std::vector<int> allEdgeIndicesRaw(1); // all indices belonging to the edge, out-of-sequence
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
                        
                        if (neighbourXPos < 0 || neighbourXPos >= haarWidth || neighbourYPos < 0 || neighbourYPos >= haarWidth)
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
            std::vector<int> branchSizes; // record length of each branch
            std::vector<std::vector<int>> branchIndicesVectors; // record indices of each branch
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

                        if (neighbourXPos < 0 || neighbourXPos >= haarWidth || neighbourYPos < 0 || neighbourYPos >= haarWidth)
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

        std::vector<double> edgePointCurvatures(edgeLength);

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
            int iStartBreakPoint = breakPointIndices[iBreakPoint] + 1;
            int iEndBreakPoint = breakPointIndices[iBreakPoint + 1];
            
            int partialEdgeLength = iEndBreakPoint - iStartBreakPoint;
            
            if (partialEdgeLength < curvatureWindowLength)
            {
                continue; // ignore short edges
            }

            // grab indices of edge points

            std::vector<int> partialEdgeIndices(allEdgeIndices.begin() + iStartBreakPoint, allEdgeIndices.begin() + iEndBreakPoint);

            // calculate pixel intensities within inner curve of edge points

            std::vector<double> partialEdgeIntensities(partialEdgeLength);

            {
                uchar *ptr_img = img.data;

                for (int iEdgePoint = 0; iEdgePoint < partialEdgeLength; iEdgePoint++)
                {
                    int edgePointIndex = partialEdgeIndices[iEdgePoint];
                    int edgePointXPos = edgePointIndex % haarWidth;
                    int edgePointYPos = (edgePointIndex - edgePointXPos) / haarWidth;

                    int offsetXPos = edgePointXPos + edgeIntensityPositionOffset * ceil2(edgePointXNormals[iEdgePoint]);
                    int offsetYPos = edgePointYPos + edgeIntensityPositionOffset * ceil2(edgePointYNormals[iEdgePoint]);

                    if (offsetXPos < 0 || offsetXPos > haarWidth || offsetYPos < 0 || offsetYPos > haarWidth)
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

            std::vector<double> partialEdgeCurvatures(edgePointCurvatures.begin() + iStartBreakPoint, edgePointCurvatures.begin() + iEndBreakPoint);
            double avgCurvature = 0;

            for (int iEdgePoint = 0; iEdgePoint < partialEdgeLength; iEdgePoint++)
            {
                avgCurvature += std::abs(partialEdgeCurvatures[iEdgePoint] / partialEdgeLength);
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

                    if (neighbourXPos < 0 || neighbourXPos >= haarWidth || neighbourYPos < 0 || neighbourYPos >= haarWidth)
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
            
            int numberOfEdgePoints = partialEdgeIndices.size();

            mEdgePropertiesAll.curvatures.push_back(avgCurvature);
            mEdgePropertiesAll.intensities.push_back(avgIntensity);
            mEdgePropertiesAll.lengths.push_back(partialEdgeLength); // add edge length to structure
            mEdgePropertiesAll.sizes.push_back(numberOfEdgePoints); // add total number of edge points of edge to structure
            mEdgePropertiesAll.indices.push_back(partialEdgeIndices);
        }
    }

    return mEdgePropertiesAll; // return structure
}

edgeProperties edgeThreshold(eyeProperties mEyeProperties, const edgeProperties& mEdgePropertiesAll, int haarWidth, double expectedCurvature)
{
    //    double edgeIntensityThreshold = mEyeProperties.v.edgeIntensity + mEyeProperties.p.edgeIntensityOffset;
    //    int maximumNumberEdges = mEyeProperties.p.edgeMaximumFitNumber;

    //    edgeProperties mEdgePropertiesNew;
    //    int numberOfEdgesTotal = mEdgePropertiesAll.sizes.size();

    //    for (int iEdge = 0; iEdge < numberOfEdgesTotal; iEdge++) // THRESHOLD: ignore bright edges
    //    {
    //        if (mEdgePropertiesAll.intensities[iEdge] < edgeIntensityThreshold)
    //        {
    //            mEdgePropertiesNew.lengths.push_back(mEdgePropertiesAll.lengths[iEdge]);
    //            mEdgePropertiesNew.sizes.push_back(mEdgePropertiesAll.sizes[iEdge]);
    //            mEdgePropertiesNew.indices.push_back(mEdgePropertiesAll.indices[iEdge]);
    //            mEdgePropertiesNew.intensities.push_back(mEdgePropertiesAll.intensities[iEdge]);
    //        }
    //    }

    int maximumNumberEdges = mEyeProperties.p.edgeMaximumFitNumber;
    edgeProperties mEdgePropertiesNew = mEdgePropertiesAll;
    int numberOfEdgesTotal = mEdgePropertiesAll.sizes.size();

    if (numberOfEdgesTotal > maximumNumberEdges) // THRESHOLD: do not handle more than a fixed number of edges
    {
        edgeProperties mEdgePropertiesFinal;

        // calculate distance between each edge point and expected pupil centre

        int numberOfEdges = mEdgePropertiesNew.sizes.size();

        double centreXPos = mEyeProperties.v.xPosPredicted - mEyeProperties.m.offsetPupilHaarXPos;
        double centreYPos = mEyeProperties.v.yPosPredicted - mEyeProperties.m.offsetPupilHaarYPos;

        std::vector<double> dRTotalsUnsorted(numberOfEdges, 0);

        for (int iEdge = 0; iEdge < numberOfEdges; iEdge++)
        {
            for (int iEdgePoint = 0, edgeSize = mEdgePropertiesNew.sizes[iEdge]; iEdgePoint < edgeSize; iEdgePoint++)
            {
                int edgePointIndex = mEdgePropertiesNew.indices[iEdge][iEdgePoint];
                int edgePointXPos = edgePointIndex % haarWidth;
                int edgePointYPos = (edgePointIndex - edgePointXPos) / haarWidth;

                double dX = centreXPos - edgePointXPos;
                double dY = centreYPos - edgePointYPos;

                double distance = sqrt(pow(dX, 2) + pow(dY, 2));
                double dR = distance - 0.75 * mEyeProperties.v.pupilRadius;
                dRTotalsUnsorted[iEdge] += dR / edgeSize; // additional weight to longer edges
            }
        }

        std::vector<double> totalWeightsUnsorted(numberOfEdges);

        for (int iEdge = 0; iEdge < numberOfEdges; iEdge++)
        {
            // intensity weight
            double intensity = mEdgePropertiesAll.intensities[iEdge];
            double intensityWeight = (20 / (1 + 0.01 * pow(0.90, -intensity + mEyeProperties.v.edgeIntensity)));
            if (intensityWeight < 0) { intensityWeight = 0; }

            // position weight
            double dR = dRTotalsUnsorted[iEdge];
            double positionWeight = (-15 / (mEyeProperties.v.pupilRadius)) * dR + 15;
            if (positionWeight < 0) { positionWeight = 0; }

            // length weight
            double length = mEdgePropertiesAll.lengths[iEdge];
            double lengthWeight;
            if (length <= mEyeProperties.v.pupilCircumference)
            {
                lengthWeight = 12 * (1 - exp(-0.0002 * mEyeProperties.v.pupilCircumference * length));
            }
            else
            {
                lengthWeight = 12 / (1 + 0.01 * pow(0.85, -length + mEyeProperties.v.pupilCircumference));
            }
            if (lengthWeight < 0) { intensityWeight = 0; }

            // curvature weight
            double dC = std::abs(mEdgePropertiesAll.curvatures[iEdge] - expectedCurvature);
            double curvatureWeight = (-7 / expectedCurvature) * dC + 7;
            if (curvatureWeight < 0) { curvatureWeight = 0; }

            totalWeightsUnsorted[iEdge] = intensityWeight + positionWeight + lengthWeight + curvatureWeight;
        }

        // Pick edges that resemble expected pupil edge the most, based on their position

        std::vector<double> totalWeightsSorted = totalWeightsUnsorted;
        std::sort(totalWeightsSorted.begin(), totalWeightsSorted.end());
        std::reverse(totalWeightsSorted.begin(), totalWeightsSorted.end());

        for (int iEdge = 0; iEdge < maximumNumberEdges; iEdge++)
        {
            for (int jEdge = 0; jEdge < numberOfEdges; jEdge++)
            {
                if (totalWeightsSorted[iEdge] == totalWeightsUnsorted[jEdge])
                {
                    mEdgePropertiesFinal.lengths.push_back(mEdgePropertiesNew.lengths[jEdge]);
                    mEdgePropertiesFinal.sizes.push_back(mEdgePropertiesNew.sizes[jEdge]);
                    mEdgePropertiesFinal.indices.push_back(mEdgePropertiesNew.indices[jEdge]);
                    mEdgePropertiesFinal.intensities.push_back(mEdgePropertiesNew.intensities[jEdge]);
                    break;
                }
            }
        }

        return mEdgePropertiesFinal;
    }

    return mEdgePropertiesNew;
}

std::vector<double> EllipseRotationTransformation(const std::vector<double>& c)
{
    double A = c[0];
    double B = c[1];
    double C = c[2];
    double D = c[3];
    double E = c[4];
    double F = c[5];
    
    double alpha = 0.5 * (atan(B / (A - C)));
    
    double AA = A * cos(alpha) * cos(alpha) + B * cos(alpha) * sin(alpha) + C * sin(alpha) * sin(alpha);
    double CC = A * sin(alpha) * sin(alpha) - B * cos(alpha) * sin(alpha) + C * cos(alpha) * cos(alpha);
    double DD = D * cos(alpha) + E * sin(alpha);
    double EE = -D * sin(alpha) + E * cos(alpha);
    double FF = F;
    
    /// semi axes
    
    double a = sqrt((-4 * FF * AA * CC + CC * DD * DD + AA * EE * EE)/(4 * AA * CC * CC));
    double b = sqrt((-4 * FF * AA * CC + CC * DD * DD + AA * EE * EE)/(4 * AA * AA * CC));
    
    /// coordinates of centre point
    
    double x = -(DD / (2 * AA)) * cos(alpha) + (EE / (2 * CC)) * sin(alpha);
    double y = -(DD / (2 * AA)) * sin(alpha) - (EE / (2 * CC)) * cos(alpha);
    
    std::vector<double> v(4);
    
    v[0] = a;
    v[1] = b;
    v[2] = x;
    v[3] = y;
    
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

    double minEigenValue = pow(10, 50); // arbitrarily large number
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
    {
        ellipseFitCoefficients[iCoefs] = (1 / sqrt(normalizationFactor)) * eigenVector(iCoefs);
    }

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
    {
        mEllipseProperties.pupilDetected = false;
    }

    double h = pow((majorAxis - minorAxis), 2) / pow((majorAxis + minorAxis), 2);
    mEllipseProperties.circumference = M_PI * (majorAxis + minorAxis) * (1 + (3 * h) / (10 + sqrt(4 - 3 * h))); // ramanujans 2nd approximation
    mEllipseProperties.fraction = minorAxis / majorAxis;
    mEllipseProperties.radius = 0.5 * (minorAxis + majorAxis);
    mEllipseProperties.xPos = ellipseParameters[2];
    mEllipseProperties.yPos = ellipseParameters[3];
    mEllipseProperties.coefficients = ellipseFitCoefficients;

    return mEllipseProperties;
}

ellipseProperties findBestEllipseFit(const edgeProperties& mEdgeProperties, int haarWidth, eyeProperties mEyeProperties, std::vector<int> cannyEdgeIndices)
{
    ellipseProperties mEllipseProperties;
    mEllipseProperties.pupilDetected = false;
    
    int totalNumberOfEdges = mEdgeProperties.sizes.size(); // total number of edges
    int totalNumberOfEdgePoints = cannyEdgeIndices.size();

    // vectors to record information for each accepted ellipse fit

    std::vector<double> allEllipseFitIntensities;
    std::vector<double> allEllipseFitRadii;
    std::vector<double> allEllipseFitCircumferences;
    std::vector<double> allEllipseFitFractions;
    std::vector<int> allEllipseFitLengths;
    std::vector<double> allEllipseFitXPos;
    std::vector<double> allEllipseFitYPos;
    std::vector<std::vector<double>> allEllipseFitCoefficients;
    
    int numberOfFits = 0; // iterator over all accepted combinations
    
    for (int combiNumberOfEdges = totalNumberOfEdges; combiNumberOfEdges >= 1; combiNumberOfEdges--) // loop through all possible edge set sizes
    {
        std::vector<bool> edgeCombination(totalNumberOfEdges);
        std::fill(edgeCombination.begin() + totalNumberOfEdges - combiNumberOfEdges, edgeCombination.end(), true);
        
        do // loop through all possible edge combinations for the current set size
        {
            std::vector<double> combiEdgeIntensities(combiNumberOfEdges);
            std::vector<int> combiEdgeLengths(combiNumberOfEdges);
            std::vector<int> combiEdgeSizes(combiNumberOfEdges);
            std::vector<std::vector<int>> combiEdgeIndices(combiNumberOfEdges);
            
            for (int iEdge = 0, jEdge = 0; iEdge < totalNumberOfEdges; ++iEdge)
            {
                if (edgeCombination[iEdge])
                {
                    combiEdgeIntensities[jEdge] = mEdgeProperties.intensities[iEdge];
                    combiEdgeLengths[jEdge] = mEdgeProperties.lengths[iEdge];
                    combiEdgeSizes[jEdge] = mEdgeProperties.sizes[iEdge];
                    combiEdgeIndices[jEdge] = mEdgeProperties.indices[iEdge];
                    jEdge++;
                }
            }
            
            // calculate lengths

            int edgeSetLength = std::accumulate(combiEdgeLengths.begin(), combiEdgeLengths.end(), 0);
            
            if (edgeSetLength < mEyeProperties.v.pupilCircumference * edgeCollectionFraction)
            {
                continue; // ignore edge collections that are too short
            }

            int edgeSetSize = std::accumulate(combiEdgeSizes.begin(), combiEdgeSizes.end(), 0);

            // concatenate index vectors

            std::vector<int> edgeIndices; // vector containing all indices for fit
            edgeIndices.reserve(edgeSetSize); // preallocate memory

            for (int iEdge = 0; iEdge < combiNumberOfEdges; iEdge++)
            {
                edgeIndices.insert(edgeIndices.end(), combiEdgeIndices[iEdge].begin(), combiEdgeIndices[iEdge].end());
            }

            // fit ellipse

            ellipseProperties mEllipsePropertiesNew = fitEllipse(edgeIndices, edgeSetSize, haarWidth);

            if (!mEllipsePropertiesNew.pupilDetected) // error
            {
                continue;
            }

            // Size and shape filters
            
            if (mEllipsePropertiesNew.circumference > mEyeProperties.p.pupilCircumferenceMax || mEllipsePropertiesNew.circumference < mEyeProperties.p.pupilCircumferenceMin) // no large or small pupils
            {
                continue;
            }

            if (mEllipsePropertiesNew.fraction < mEyeProperties.p.pupilFractMin) // no extreme deviations from circular shape
            {
                continue;
            }

            if (std::abs(mEllipsePropertiesNew.circumference - mEyeProperties.v.pupilCircumference) > mEyeProperties.v.thresholdCircumferenceChange) // no large pupil size changes
            {
                continue;
            }

            if (std::abs(mEllipsePropertiesNew.fraction - mEyeProperties.v.pupilFraction) > mEyeProperties.v.thresholdFractionChange) // no large pupil size changes
            {
                continue;
            }

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
                    int edgePointIndex = cannyEdgeIndices[iEdgePoint];

                    double x = edgePointIndex % haarWidth;
                    double y = (edgePointIndex - x) / haarWidth;

                    double fitError = std::abs(A * x * x + B * x * y + C * y * y + D * x + E * y + F);

                    if (fitError < mEyeProperties.p.ellipseFitErrorMaximum)
                    {
                        fittedEdgePoints.push_back(edgePointIndex);
                    }
                }

                fitSetLength = fittedEdgePoints.size();

                if (fitSetLength < mEyeProperties.v.pupilCircumference * minimumFitFraction)
                {
                    continue; // Threshold: ignore edge collections that are too short
                }
            }

            // do refinement if ellipse fit overlaps with many new edge points

            if ((fitSetLength - edgeSetLength) > mEyeProperties.v.pupilCircumference * minimumRefinementFraction)
            {
                mEllipsePropertiesNew = fitEllipse(fittedEdgePoints, fitSetLength, haarWidth);
            }

            // save parameters of accepted fit

            allEllipseFitIntensities.push_back(calculateMean(combiEdgeIntensities));
            allEllipseFitCircumferences.push_back(mEllipsePropertiesNew.circumference);
            allEllipseFitFractions.push_back(mEllipsePropertiesNew.fraction);
            allEllipseFitRadii.push_back(mEllipsePropertiesNew.radius);
            allEllipseFitLengths.push_back(fitSetLength);
            allEllipseFitXPos.push_back(mEllipsePropertiesNew.xPos);
            allEllipseFitYPos.push_back(mEllipsePropertiesNew.yPos);
            allEllipseFitCoefficients.push_back(mEllipsePropertiesNew.coefficients);
            
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
            double circumferenceChange = (std::abs(allEllipseFitCircumferences[iFit] - mEyeProperties.v.pupilCircumference));

            if (circumferenceChange < 0.5 * mEyeProperties.v.thresholdCircumferenceChange)
            {
                circumferenceChange = 0.5 * mEyeProperties.v.thresholdCircumferenceChange;
            }

            double fractionChange = (std::abs(allEllipseFitFractions[iFit] - mEyeProperties.v.pupilFraction));

            if (fractionChange < 0.5 * mEyeProperties.v.thresholdFractionChange)
            {
                fractionChange = 0.5 * mEyeProperties.v.thresholdFractionChange;
            }

            double normalizationConstant = mEyeProperties.v.thresholdCircumferenceChange / mEyeProperties.v.thresholdFractionChange;
            fractionChange = normalizationConstant * fractionChange;

            featureChange[iFit] = (circumferenceChange + fractionChange) / sqrt((double) allEllipseFitLengths[iFit]);
        }

        int acceptedFitIndex = std::distance(featureChange.begin(), std::min_element(featureChange.begin(), featureChange.end()));

        mEllipseProperties.intensity = allEllipseFitIntensities[acceptedFitIndex];
        mEllipseProperties.circumference = allEllipseFitCircumferences[acceptedFitIndex];
        mEllipseProperties.fraction = allEllipseFitFractions[acceptedFitIndex];
        mEllipseProperties.radius = allEllipseFitRadii[acceptedFitIndex];
        mEllipseProperties.xPos = allEllipseFitXPos[acceptedFitIndex];
        mEllipseProperties.yPos = allEllipseFitYPos[acceptedFitIndex];
        mEllipseProperties.coefficients = allEllipseFitCoefficients[acceptedFitIndex];
        mEllipseProperties.pupilDetected = true;
    }
    
    return mEllipseProperties;
}

eyeProperties pupilDetector(const cv::Mat& imageOriginalBGR, eyeProperties mEyeProperties)
{
    // Define some variables
    
    eyeProperties mEyePropertiesNew = mEyeProperties; // new properties for new frame
    edgeProperties mEdgePropertiesNew;
    ellipseProperties mEllipseProperties;

    mEllipseProperties.pupilDetected = false;
    mEyePropertiesNew.m.errorDetected = false;
    mEyePropertiesNew.m.image = imageOriginalBGR;

    // Check limits

    if (mEyeProperties.v.thresholdCircumferenceChange > mEyeProperties.p.pupilCircumferenceMax)
    {
        mEyeProperties.v.thresholdCircumferenceChange = mEyeProperties.p.pupilCircumferenceMax;
    }
    else if (mEyeProperties.v.thresholdCircumferenceChange < mEyeProperties.p.thresholdCircumferenceChangeMin)
    {
        mEyeProperties.v.thresholdCircumferenceChange = mEyeProperties.p.thresholdCircumferenceChangeMin;
    }

    if (mEyeProperties.v.thresholdFractionChange > 1.0)
    {
        mEyeProperties.v.thresholdFractionChange = 1.0;
    }
    else if (mEyeProperties.v.thresholdFractionChange < mEyeProperties.p.thresholdFractionChangeMin)
    {
        mEyeProperties.v.thresholdFractionChange = mEyeProperties.p.thresholdFractionChangeMin;
    }

    // Define search area

    int imageWdth = imageOriginalBGR.cols;
    int imageHght = imageOriginalBGR.rows;

    int searchStartX = round(mEyeProperties.v.xPosPredicted - mEyeProperties.v.searchRadius);
    int searchStartY = round(mEyeProperties.v.yPosPredicted - mEyeProperties.v.searchRadius);
    
    if (searchStartX < 0)
    {
        searchStartX = 0;
    }
    
    if (searchStartY < 0)
    {
        searchStartY = 0;
    }
    
    int searchEndX = round(mEyeProperties.v.xPosPredicted + mEyeProperties.v.searchRadius);
    int searchEndY = round(mEyeProperties.v.yPosPredicted + mEyeProperties.v.searchRadius);
    
    if (searchEndX >= imageWdth)
    {
        searchEndX = imageWdth - 1;
    }
    
    if (searchEndY >= imageHght)
    {
        searchEndY = imageHght - 1;
    }
    
    int searchWdth = searchEndX - searchStartX + 1;
    int searchHght = searchEndY - searchStartY + 1;
    
    int pupilHaarWdth = pupilHaarFraction * (round(mEyeProperties.v.pupilCircumference / M_PI));

    int offsetPupilHaarXPos = 0;
    int offsetPupilHaarYPos = 0;
    int offsetPupilHaarWdth = 0;

    if (pupilHaarWdth > 0 && searchWdth > pupilHaarWdth && searchHght > pupilHaarWdth)
    {
        // Convert to grayscale

        cv::Mat imageOriginalGray;
        cv::cvtColor(imageOriginalBGR, imageOriginalGray, cv::COLOR_BGR2GRAY);

        // Initial approximate detection

        std::vector<unsigned int> integralImage = calculateIntImg(imageOriginalGray, imageWdth, searchStartX, searchStartY, searchWdth, searchHght);

        haarProperties glintHaarProperties = GlintDetector(imageOriginalGray, imageWdth, searchStartX, searchStartY, searchWdth, searchHght, 2 * mEyeProperties.p.glintRadius, mEyeProperties.p.glintRadius);

        int glintHaarXPos = searchStartX + glintHaarProperties.x_pos;
        int glintHaarYPos = searchStartY + glintHaarProperties.y_pos;

        haarProperties pupilHaarProperties = PupilHaarDetector(integralImage, searchWdth, searchHght, pupilHaarWdth, glintHaarXPos, glintHaarYPos, mEyeProperties.p.glintRadius);

        int pupilHaarXPos = searchStartX + pupilHaarProperties.x_pos;
        int pupilHaarYPos = searchStartY + pupilHaarProperties.y_pos;

        offsetPupilHaarXPos = pupilHaarXPos - round(pupilHaarWdth * mEyeProperties.p.pupilOffset);
        offsetPupilHaarYPos = pupilHaarYPos - round(pupilHaarWdth * mEyeProperties.p.pupilOffset);
        offsetPupilHaarWdth = pupilHaarWdth + round(pupilHaarWdth * mEyeProperties.p.pupilOffset * 2);

        // Check limits

        if (offsetPupilHaarWdth >= imageHght)
        {
            offsetPupilHaarWdth = imageHght - 1;
        }

        if (offsetPupilHaarXPos < 0)
        {
            offsetPupilHaarXPos = 0;
        }

        if (offsetPupilHaarYPos < 0)
        {
            offsetPupilHaarYPos = 0;
        }

        if (offsetPupilHaarXPos + offsetPupilHaarWdth >= imageWdth)
        {
            offsetPupilHaarXPos = imageWdth - offsetPupilHaarWdth - 1;
        }

        if (offsetPupilHaarYPos + offsetPupilHaarWdth >= imageHght)
        {
            offsetPupilHaarYPos = imageHght - offsetPupilHaarWdth - 1;
        }

        // Crop image to outer pupil Haar region

        cv::Rect pupilROI(offsetPupilHaarXPos, offsetPupilHaarYPos, offsetPupilHaarWdth, offsetPupilHaarWdth);
        cv::Mat imagePupilBGR = imageOriginalBGR(pupilROI);

        // Convert back to grayscale

        cv::Mat imagePupilGray;
        cv::cvtColor(imagePupilBGR, imagePupilGray, cv::COLOR_BGR2GRAY);

        // Canny edge detection

        cv::Mat imagePupilGrayBlurred;

        if (mEyeProperties.p.cannyBlurLevel > 0)
        {
            cv::GaussianBlur(imagePupilGray, imagePupilGrayBlurred, cv::Size(mEyeProperties.p.cannyBlurLevel, mEyeProperties.p.cannyBlurLevel), 0, 0);
        }
        else
        {
            imagePupilGrayBlurred = imagePupilGray;
        }

        cv::Mat imageCannyEdges;
        cv::Canny(imagePupilGrayBlurred, imageCannyEdges, mEyeProperties.p.cannyLowerLimit, mEyeProperties.p.cannyUpperLimit, mEyeProperties.p.cannyKernelSize);

        std::vector<char> cannyEdgesRaw = cannyConversion(imageCannyEdges, offsetPupilHaarWdth);
        std::vector<char> cannyEdges = MorphOpen(cannyEdgesRaw, offsetPupilHaarWdth);
        std::vector<int> cannyEdgeIndices = findCannyIndices(cannyEdges);

        // Edge thresholding

        double curvatureLowerLimit;

        {
            double A =  5796;
            double B = -0.7855;
            double C = -0.7259;
            double D =  0.3801;
            double E = -1.408;

            curvatureLowerLimit = A * exp(B * pow(mEyeProperties.v.pupilCircumference, D) + C * pow(mEyeProperties.v.pupilFraction, E)) - mEyeProperties.p.curvatureOffset;
        }

        double curvatureUpperLimit;

        {
            double A =  98.31;
            double B =  9.95;
            double C = -2.106;
            double D = -0.412;
            double E =  0.4122;

            curvatureUpperLimit = A * exp(B * pow(mEyeProperties.v.pupilCircumference, D) + C * pow(mEyeProperties.v.pupilFraction, E)) + mEyeProperties.p.curvatureOffset;
        }

        edgeProperties mEdgeProperties = EdgeFilter(imagePupilGray, cannyEdges, offsetPupilHaarWdth, curvatureLowerLimit, curvatureUpperLimit);

        if (mEdgeProperties.indices.size() > 0) // THRESHOLD: ignore empty edge collections
        {
            double expectedCurvature = 0.5 * (curvatureUpperLimit + curvatureLowerLimit);
            mEdgePropertiesNew = edgeThreshold(mEyeProperties, mEdgeProperties, offsetPupilHaarWdth, expectedCurvature);

            // Ellipse fitting

            mEllipseProperties = findBestEllipseFit(mEdgePropertiesNew, offsetPupilHaarWdth, mEyeProperties, cannyEdgeIndices);
        }

        // Save parameters

        mEyePropertiesNew.v.pupilCircumferenceExact = mEllipseProperties.circumference;
        mEyePropertiesNew.v.pupilFractionExact = mEllipseProperties.fraction;

        mEyePropertiesNew.v.pupilDetected = mEllipseProperties.pupilDetected;

        // For draw functions

        mEyePropertiesNew.m.offsetPupilHaarXPos = offsetPupilHaarXPos;
        mEyePropertiesNew.m.offsetPupilHaarYPos = offsetPupilHaarYPos;
        mEyePropertiesNew.m.offsetPupilHaarWdth = offsetPupilHaarWdth;

        mEyePropertiesNew.m.pupilHaarXPos = pupilHaarXPos;
        mEyePropertiesNew.m.pupilHaarYPos = pupilHaarYPos;
        mEyePropertiesNew.m.pupilHaarWdth = pupilHaarWdth;

        mEyePropertiesNew.m.glintHaarXPos = glintHaarXPos;
        mEyePropertiesNew.m.glintHaarYPos = glintHaarYPos;
        mEyePropertiesNew.m.glintHaarWidth = 2 * mEyeProperties.p.glintRadius;

        mEyePropertiesNew.m.cannyEdges = cannyEdgesRaw;
        mEyePropertiesNew.m.edgeIndices = mEdgeProperties.indices;
        mEyePropertiesNew.m.edgeIndicesNew = mEdgePropertiesNew.indices;
        mEyePropertiesNew.m.ellipseCoefficients = mEllipseProperties.coefficients;
    }
    else
    {
        mEyePropertiesNew.m.errorDetected = true;
    }

    // For running averages
    
    if (!mEllipseProperties.pupilDetected)
    {
        mEyePropertiesNew.v.pupilFraction = mEyeProperties.v.pupilFraction + mEyeProperties.p.alphaPupil * (mEyeProperties.p.pupilFractIni - mEyeProperties.v.pupilFraction);
        mEyePropertiesNew.v.pupilCircumference = mEyeProperties.v.pupilCircumference + mEyeProperties.p.alphaPupil * (mEyeProperties.p.pupilCircfIni - mEyeProperties.v.pupilCircumference);
        mEyePropertiesNew.v.pupilRadius = mEyeProperties.v.pupilCircumference / (2 * M_PI);

        mEyePropertiesNew.v.momentumFraction = mEyeProperties.v.momentumFraction * mEyeProperties.p.alphaMomentum;
        mEyePropertiesNew.v.momentumCircumference = mEyeProperties.v.momentumCircumference * mEyeProperties.p.alphaMomentum;
        mEyePropertiesNew.v.momentumRadius = mEyeProperties.v.momentumRadius * mEyeProperties.p.alphaMomentum;
        
        mEyePropertiesNew.v.edgeIntensity = mEyeProperties.v.edgeIntensity + mEyeProperties.p.alphaPupil * (mEyeProperties.p.edgeIntensityIni - mEyeProperties.v.edgeIntensity);

        mEyePropertiesNew.v.xVelocity = mEyeProperties.v.xVelocity * mEyeProperties.p.alphaVelocity;
        mEyePropertiesNew.v.yVelocity = mEyeProperties.v.yVelocity * mEyeProperties.p.alphaVelocity;

        mEyePropertiesNew.v.xPosPredicted = mEyeProperties.v.xPosPredicted + mEyeProperties.p.alphaPosition * (offsetPupilHaarXPos + 0.5 * offsetPupilHaarWdth - mEyeProperties.v.xPosPredicted) + mEyeProperties.v.xVelocity;
        mEyePropertiesNew.v.yPosPredicted = mEyeProperties.v.yPosPredicted + mEyeProperties.p.alphaPosition * (offsetPupilHaarYPos + 0.5 * offsetPupilHaarWdth - mEyeProperties.v.yPosPredicted) + mEyeProperties.v.yVelocity;

        mEyePropertiesNew.v.searchRadius = mEyeProperties.v.searchRadius * (2 - mEyeProperties.p.alphaGeneral);
        mEyePropertiesNew.v.thresholdCircumferenceChange = mEyeProperties.v.thresholdCircumferenceChange * (2 - mEyeProperties.p.alphaGeneral);
        mEyePropertiesNew.v.thresholdFractionChange = mEyeProperties.v.thresholdFractionChange * (2 - mEyeProperties.p.alphaGeneral);
    }
    else
    {
        mEyePropertiesNew.v.xPos = mEllipseProperties.xPos + offsetPupilHaarXPos;
        mEyePropertiesNew.v.yPos = mEllipseProperties.yPos + offsetPupilHaarYPos;

        mEyePropertiesNew.v.pupilFraction = mEyeProperties.v.pupilFraction + mEyeProperties.p.alphaPupil * (mEllipseProperties.fraction - mEyeProperties.v.pupilFraction) + mEyeProperties.v.momentumFraction;
        mEyePropertiesNew.v.pupilCircumference = mEyeProperties.v.pupilCircumference + mEyeProperties.p.alphaPupil * (mEllipseProperties.circumference - mEyeProperties.v.pupilCircumference) + mEyeProperties.v.momentumCircumference;
        mEyePropertiesNew.v.pupilRadius = mEyeProperties.v.pupilRadius + mEyeProperties.p.alphaPupil * (mEllipseProperties.radius - mEyeProperties.v.pupilRadius) + mEyeProperties.v.momentumRadius;

        mEyePropertiesNew.v.momentumFraction = (mEyeProperties.v.momentumFraction + (mEyePropertiesNew.v.pupilFraction - mEyeProperties.v.pupilFraction)) * mEyeProperties.p.alphaMomentum;
        mEyePropertiesNew.v.momentumCircumference = (mEyeProperties.v.momentumCircumference + (mEyePropertiesNew.v.pupilCircumference - mEyeProperties.v.pupilCircumference)) * mEyeProperties.p.alphaMomentum;
        mEyePropertiesNew.v.momentumRadius = (mEyeProperties.v.momentumRadius + (mEyePropertiesNew.v.pupilRadius - mEyeProperties.v.pupilRadius)) * mEyeProperties.p.alphaMomentum;

        mEyePropertiesNew.v.edgeIntensity = mEyeProperties.v.edgeIntensity + mEyeProperties.p.alphaPupil * (mEllipseProperties.intensity - mEyeProperties.v.edgeIntensity);
        
        mEyePropertiesNew.v.xVelocity = (mEyeProperties.v.xVelocity + (mEyePropertiesNew.v.xPosPredicted - mEyeProperties.v.xPosPredicted)) * mEyeProperties.p.alphaVelocity;
        mEyePropertiesNew.v.yVelocity = (mEyeProperties.v.yVelocity + (mEyePropertiesNew.v.yPosPredicted - mEyeProperties.v.yPosPredicted)) * mEyeProperties.p.alphaVelocity;

        mEyePropertiesNew.v.xPosPredicted = mEyeProperties.v.xPosPredicted + mEyeProperties.p.alphaPosition * (mEyePropertiesNew.v.xPos - mEyeProperties.v.xPosPredicted) + mEyeProperties.v.xVelocity;
        mEyePropertiesNew.v.yPosPredicted = mEyeProperties.v.yPosPredicted + mEyeProperties.p.alphaPosition * (mEyePropertiesNew.v.yPos - mEyeProperties.v.yPosPredicted) + mEyeProperties.v.yVelocity;
        
        mEyePropertiesNew.v.searchRadius = mEyeProperties.v.searchRadius * mEyeProperties.p.alphaGeneral;
        mEyePropertiesNew.v.thresholdCircumferenceChange = mEyeProperties.v.thresholdCircumferenceChange * mEyeProperties.p.alphaGeneral;
        mEyePropertiesNew.v.thresholdFractionChange = mEyeProperties.v.thresholdFractionChange * mEyeProperties.p.alphaGeneral;
    }
    
    if (mEyePropertiesNew.v.searchRadius > imageWdth)
    {
        mEyePropertiesNew.v.searchRadius = imageWdth;
    }
    else if (mEyePropertiesNew.v.searchRadius < (0.5 * offsetPupilHaarWdth))
    {
        mEyePropertiesNew.v.searchRadius = ceil(0.5 * offsetPupilHaarWdth);
    }

    if (mEyePropertiesNew.v.thresholdCircumferenceChange > mEyePropertiesNew.p.pupilCircumferenceMax)
    {
        mEyePropertiesNew.v.thresholdCircumferenceChange = mEyePropertiesNew.p.pupilCircumferenceMax;
    }
    else if (mEyePropertiesNew.v.thresholdCircumferenceChange < mEyePropertiesNew.p.thresholdCircumferenceChangeMin)
    {
        mEyePropertiesNew.v.thresholdCircumferenceChange = mEyePropertiesNew.p.thresholdCircumferenceChangeMin;
    }

    if (mEyePropertiesNew.v.thresholdFractionChange > 1.0)
    {
        mEyePropertiesNew.v.thresholdFractionChange = 1.0;
    }
    else if (mEyePropertiesNew.v.thresholdFractionChange < mEyePropertiesNew.p.thresholdFractionChangeMin)
    {
        mEyePropertiesNew.v.thresholdFractionChange = mEyePropertiesNew.p.thresholdFractionChangeMin;
    }

    return mEyePropertiesNew;
}

bool flashDetection(const cv::Mat& imgBGR, double threshold)
{
    cv::Mat img;
    cv::cvtColor(imgBGR, img, cv::COLOR_BGR2GRAY);

    bool flashDetected = false;
    unsigned long long intensityTotal = 0;
    uchar *ptr = img.data;
    int imgSize = img.cols * img.rows;

    for (int iPixel = 0; iPixel < imgSize; iPixel++)
    {
        intensityTotal += ptr[iPixel];
    }

    double avgIntensity = intensityTotal / (double) imgSize;

    if (avgIntensity > (255 * threshold))
    {
        flashDetected = true; // flash detected
    }

    return flashDetected;
}
