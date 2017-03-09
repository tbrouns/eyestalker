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

int find2(std::vector<std::vector<int>>& v, int target)
{
    int size = v.size();
    for (int i = 0; i < size; i++)
    {
        auto itr = find(v[i].begin(), v[i].end(), target);
        if (v[i].end() != itr) { return i; }
    }
    return size;
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

inline double gaussian(double x, double sigma)
{
    return exp(-pow(x/sigma,2));
}

double calculateScoreTotal(const detectionVariables& mDetectionVariables, std::vector<double>& featureValues, bool USE_LENGTH, bool USE_CERTAINTY)
{   
    static const double weightCircumference = 0.24;
    static const double weightIntensity     = 0.66;
    static const double weightCurvature     = 0.42;
    static const double weightRadius        = 0.65;
    static const double weightRadiusVar     = 0.21;
    static const double weightGradient      = 0.24;
    static const double weightBeta          = 0.49;
    
    static const double sigmaRadius        = 0.0655810;
    static const double sigmaCircumference = 0.7269800;
    static const double sigmaCurvature     = 3.4627000;
    static const double sigmaIntensity     = 11.532900;
    static const double sigmaGradient      = 5.3539000;
    static const double sigmaRadiusVar     = 0.0039498;
    
    for (int i = 0, vSize = featureValues.size(); i < vSize; i++) // check for NaNs or Infs
    {
        double val = featureValues[i];
        if (!std::isfinite(val)) { featureValues[i] = 0; }
        else if (val < 0)        { featureValues[i] = std::abs(val); }
    }
    
    // Do score calculation
    
    double featureValueRadius        = featureValues[0];
    double featureValueCircumference = featureValues[1];
    double featureValueCurvature     = featureValues[2];
    double featureValueIntensity     = featureValues[3];
    double featureValueGradient      = featureValues[4];
    double featureValueRadiusVar     = featureValues[5];
    
    double certaintyFactorPosition = mDetectionVariables.certaintyPosition;
    double certaintyFactorFeatures = mDetectionVariables.certaintyFeatures;
    
    double factorRadius        = certaintyFactorPosition  * weightRadius;
    double factorCircumference = certaintyFactorFeatures  * weightCircumference;
    
    double certaintyFactorIntensity;
    double certaintyFactorCurvature;
    double certaintyFactorGradient;
    
    if (!USE_LENGTH) { factorCircumference = 0; }
    
    if (!USE_CERTAINTY) // true when comparing feature values with predicted values
    {
        certaintyFactorIntensity = 1.0;
        certaintyFactorCurvature = 1.0;
        certaintyFactorGradient  = certaintyFactorPosition;
    }
    else
    {
        certaintyFactorIntensity = certaintyFactorFeatures;
        certaintyFactorCurvature = certaintyFactorFeatures;
        certaintyFactorGradient  = certaintyFactorPosition * certaintyFactorFeatures;
    }
    
    double factorIntensity = certaintyFactorIntensity * weightIntensity;
    double factorGradient  = certaintyFactorGradient  * weightGradient;
    
    // Importance of curvature and radial variance should be dependent on length of edge
    
    double factorLength = weightBeta * featureValueCircumference + (1 - weightBeta);
    
    double factorRadiusVar = certaintyFactorPosition  * weightRadiusVar * factorLength;
    double factorCurvature = certaintyFactorCurvature * weightCurvature * factorLength;
    
    // Calculate scores
    
    double scoreRadius        = factorRadius        * gaussian(featureValueRadius,        sigmaRadius       );
    double scoreCircumference = factorCircumference * gaussian(featureValueCircumference, sigmaCircumference);
    double scoreCurvature     = factorCurvature     * gaussian(featureValueCurvature,     sigmaCurvature    );
    double scoreIntensity     = factorIntensity     * gaussian(featureValueIntensity,     sigmaIntensity    );
    double scoreGradient      = factorGradient      * gaussian(featureValueGradient,      sigmaGradient     );
    double scoreRadiusVar     = factorRadiusVar     * gaussian(featureValueRadiusVar,     sigmaRadiusVar    );
    
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

std::vector<int> getEdgeIndices(const std::vector<int>& binaryImageVector, int tag)
{
    int AOIArea = binaryImageVector.size();
    std::vector<int> cannyEdgeIndices;
    for (int iEdgePoint = 0; iEdgePoint < AOIArea; iEdgePoint++)
    { if (binaryImageVector[iEdgePoint] == tag) { cannyEdgeIndices.push_back(iEdgePoint); }}
    return cannyEdgeIndices;
}

std::vector<int> sharpenEdges_1(std::vector<int>& binaryImageVector, std::vector<int>& edgePointIndicesOld, AOIProperties mAOI)
{
    std::vector<int> edgePointIndices;
    
    // First morphological operation
    int numEdgePoints = edgePointIndicesOld.size();
    
    std::vector<int> dX = {  0,  1, -1,  1,  0, -1,  1, -1};
    std::vector<int> dY = { -1,  1,  0, -1,  1, -1,  0,  1};
    
    for (int iEdgePoint = 0; iEdgePoint < numEdgePoints; iEdgePoint++)
    {
        int iCentre = edgePointIndicesOld[iEdgePoint];
        int xCentre = iCentre % mAOI.wdth;
        int yCentre = (iCentre - xCentre) / mAOI.wdth;
        
        for (int m = 0; m < 4; m++)
        {
            int numFilledPixels = 0;
            // check combination of two neighbouring pixels in 4-connected environment
            for (int n = 0; n < 2; n++) // loop through two connected neighbouring pixels
            {
                int q = 2 * (m + n) % 8;
                int xNeighbour = xCentre + dX[q];
                int yNeighbour = yCentre + dY[q];
                if (xNeighbour < 0 || xNeighbour >= mAOI.wdth ||yNeighbour < 0 || yNeighbour >= mAOI.hght) { continue; } // neighbour is out-of-bounds
                int iNeighbour = mAOI.wdth * yNeighbour + xNeighbour;
                if (binaryImageVector[iNeighbour] == 1) { numFilledPixels++; } // check if neighbour is filled
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
                    break;
                }
            }
            
            edgePointIndices.push_back(iCentre); // add to vector if not removed
        }
    }
    
    return edgePointIndices;
}

std::vector<int> sharpenEdges_2(std::vector<int>& binaryImageVector, std::vector<int>& edgePointIndicesOld, AOIProperties mAOI)
{
    std::vector<int> edgePointIndices;
    
    // Second morphological operation
    
    int numEdgePoints = edgePointIndicesOld.size();
    
    std::vector<int> dX = {  1,  1,  0, -1, -1, -1,  0,  1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};
    
    for (int iEdgePoint = 0; iEdgePoint < numEdgePoints; iEdgePoint++)
    {
        int iCentre = edgePointIndicesOld[iEdgePoint];
        int xCentre = iCentre % mAOI.wdth;
        int yCentre = (iCentre - xCentre) / mAOI.wdth;
        std::vector<int> filledPixels;
        for (int m = 0; m < 8; m++)
        {
            int xNeighbour = xCentre + dX[m];
            int yNeighbour = yCentre + dY[m];
            if (xNeighbour < 0 || xNeighbour >= mAOI.wdth || yNeighbour < 0 || yNeighbour >= mAOI.hght) { continue; } // neighbour is out-of-bounds
            int iNeighbour = mAOI.wdth * yNeighbour + xNeighbour;
            if (binaryImageVector[iNeighbour] == 1) { filledPixels.push_back(m); } // check if neighbour is filled
        }
        
        int numFilledPixels = filledPixels.size();
        if (numFilledPixels == 2 || numFilledPixels == 3)
        {
            int k = 0;
            for (int i = 0; i < numFilledPixels; i++)
            {
                int j = (i + 1) % numFilledPixels;
                if (filledPixels[i] != filledPixels[j]) { k++; }
            }
            if (k == 1)
            {
                binaryImageVector[iCentre] = -1;
                continue;
            }
        }
        
        edgePointIndices.push_back(iCentre); // add to vector if not removed
    }
    
    return edgePointIndices;
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

std::vector<int> findEdges(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, std::vector<int>& cannyEdgeVector, AOIProperties mAOI)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};
    
    double radiusMax = mDetectionVariables.thresholdChangePosition + mDetectionVariables.predictedCircumference * (1 + mDetectionVariables.thresholdChangeCircumference) / (2 * M_PI);
    double radiusUpperLimit = mDetectionParameters.circumferenceMax / (2 * M_PI);
    if (radiusMax > radiusUpperLimit) { radiusMax = radiusUpperLimit; }
    
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
                if (tag == 1 || tag == 7) { EDGE_FOUND = true; }
                if (tag == 1) { startIndices.push_back(centreIndex); }
                
                if (EDGE_FOUND)
                {
                    if (R > radiusMax) { STOP_SEARCH = true; }
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
        
        double xTangent = -calculateMean(edgeXTangents); // reverse direction
        double yTangent = -calculateMean(edgeYTangents);
        
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
        
        if (neighbourXPos >= 0 && neighbourXPos < mAOI.wdth && neighbourYPos >= 0 && neighbourYPos < mAOI.hght) // check if neighbour is out-of-bounds
        {
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
    }
    
    return startIndex;
}

std::vector<int> findEdgePoints(const detectionParameters& mDetectionParameters, std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int startIndex)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};
    
    cannyEdgeVector[startIndex] = 2; // tag pixel
    
    std::vector<int> edgePointsOld = {startIndex};
    std::vector<int> edgePointsAll;
    std::vector<int> edgeTerminals;
    
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
                else if (neighbourTag == 2 || neighbourTag == 3) { nConnections++; }
            }
            
            if (nConnections == 1) // start or end of edge
            {
                edgeTerminals.push_back(centreIndex);
            }
        }
        
        numEdgePoints = edgePointsNew.size();
        
        if (numEdgePoints == 0) // check terminals once all other points have been found
        {
            int numTerminals = edgeTerminals.size();
            for (int iEdgePoint = 0; iEdgePoint < numTerminals; iEdgePoint++)
            {
                int terminalIndex = edgeTerminals[iEdgePoint];
                int edgePointNew = connectEdges(mDetectionParameters, cannyEdgeVector, mAOI, terminalIndex); // connect possible edge terminals
                
                if (terminalIndex != edgePointNew)
                {
                    cannyEdgeVector[edgePointNew] = 2; // tag newly added point
                    edgePointsNew.push_back(edgePointNew);
                }
            }
            
            numEdgePoints = edgePointsNew.size();
        }
        
        edgePointsAll.insert(edgePointsAll.begin(), edgePointsNew.begin(), edgePointsNew.end());
        edgePointsOld = edgePointsNew;
        edgePointsNew.clear();
        
    } while (numEdgePoints > 0);
    
    return edgePointsAll;
}

std::vector<vertexProperties> findGraphVertices(std::vector<int>& cannyEdgeVector, AOIProperties mAOI, int startIndex)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};
    
    cannyEdgeVector[startIndex] = 3; // tag pixel
    std::vector<int> edgePointsOld = {startIndex};
    int numEdgePoints = 0;
    std::vector<vertexProperties> verticesAll;
    std::vector<int> vertexPointsAllRaw;
    
    do
    {
        std::vector<int> edgePointsNew;
        numEdgePoints = edgePointsOld.size();
        for (int iEdgePoint = 0; iEdgePoint < numEdgePoints; iEdgePoint++) // loop through all newly added unchecked edge points
        {
            int centreIndex  = edgePointsOld[iEdgePoint]; // index of current edge point
            int centreXPos   = centreIndex % mAOI.wdth;
            int centreYPos   = (centreIndex - centreXPos) / mAOI.wdth;
            int nConnections = 0;
            
            for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
            {
                int neighbourXPos = centreXPos + dX[m];
                int neighbourYPos = centreYPos + dY[m];
                if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds
                int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;
                int neighbourTag = cannyEdgeVector[neighbourIndex];
                if (neighbourTag > 1)
                {
                    nConnections++;
                    if (neighbourTag == 2)
                    {
                        edgePointsNew.push_back(neighbourIndex);
                        cannyEdgeVector[neighbourIndex] = 3;
                    }
                }
            }
            
            // Check if a vertex was found
            if (nConnections == 1) // external vertex
            {
                cannyEdgeVector[centreIndex] = 5;
                vertexProperties vertexNew;
                vertexNew.pointIndices.push_back(centreIndex);
                vertexNew.tag = 1;
                verticesAll.push_back(vertexNew);
            }
            else if (nConnections >= 3) // internal vertex
            {
                cannyEdgeVector[centreIndex] = 4; // temporary tag
                vertexPointsAllRaw.push_back(centreIndex);
            }
        }
        edgePointsOld = edgePointsNew;
        numEdgePoints  = edgePointsOld.size();
        edgePointsNew.clear();
    } while (numEdgePoints > 0);
    
    // Create internal vertices
    
    for (int iVertex = 0, numVertices = vertexPointsAllRaw.size(); iVertex < numVertices; iVertex++)
    {
        int vertexPointCentre = vertexPointsAllRaw[iVertex];
        if (cannyEdgeVector[vertexPointCentre] == 4)
        {
            std::vector<int> vertexPointsOld = {vertexPointCentre};
            std::vector<int> vertexPointsAll = vertexPointsOld;
            cannyEdgeVector[vertexPointCentre] = 5;
            int numEdgePoints = 1;
            do
            {
                std::vector<int> vertexPointsNew;
                for (int iEdgePoint = 0; iEdgePoint < numEdgePoints; iEdgePoint++)
                {
                    int vertexPoint = vertexPointsOld[iEdgePoint];
                    int vertexXPos  = vertexPoint % mAOI.wdth;
                    int vertexYPos  = (vertexPoint - vertexXPos) / mAOI.wdth;
                    for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
                    {
                        int neighbourXPos = vertexXPos + dX[m];
                        int neighbourYPos = vertexYPos + dY[m];
                        if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds
                        int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;
                        int neighbourTag = cannyEdgeVector[neighbourIndex];
                        if (neighbourTag == 4)
                        {
                            cannyEdgeVector[neighbourIndex] = 5;
                            vertexPointsNew.push_back(neighbourIndex);
                        }
                    }
                }
                
                vertexPointsAll.insert(vertexPointsAll.end(), vertexPointsNew.begin(), vertexPointsNew.end());
                vertexPointsOld = vertexPointsNew;
                numEdgePoints = vertexPointsOld.size();
                vertexPointsNew.clear();
            } while (numEdgePoints > 0);
            
            // Add vertex
            vertexProperties vertexNew;
            vertexNew.pointIndices = vertexPointsAll;
            vertexNew.tag = 2;
            verticesAll.push_back(vertexNew);
        }
    }
    
    // Create vertex if full cyclic edge
    
    if (verticesAll.size() == 0)
    {
        vertexProperties vertexNew;
        vertexNew.pointIndices.push_back(startIndex);
        vertexNew.tag = 2;
        verticesAll.push_back(vertexNew);
    }
    
    return verticesAll;
}

std::vector<branchProperties> findGraphBranches(const detectionParameters& mDetectionParameters, std::vector<vertexProperties>& verticesAll, std::vector<int>& cannyEdgeVector, AOIProperties mAOI)
{
    std::vector<int> dX = { -1, -1,  0,  1,  1,  1,  0, -1};
    std::vector<int> dY = {  0, -1, -1, -1,  0,  1,  1,  1};
    
    std::vector<branchProperties> branchesAll;
    
    int numVertices = verticesAll.size();
    std::vector<std::vector<int>> vertexPointIndices(numVertices);
    for (int iVertex = 0; iVertex < numVertices; iVertex++)
    { vertexPointIndices[iVertex] = verticesAll[iVertex].pointIndices; }
    
    for (int iVertex = 0; iVertex < numVertices; iVertex++)
    {
        // Find all branches that are connected to vertex
        int numVertexPoints = verticesAll[iVertex].pointIndices.size();
        std::vector<int> connectedPoints;
        for (int iEdgePoint = 0; iEdgePoint < numVertexPoints; iEdgePoint++)
        {
            int centreIndex = verticesAll[iVertex].pointIndices[iEdgePoint];
            int centreXPos = centreIndex % mAOI.wdth;
            int centreYPos = (centreIndex - centreXPos) / mAOI.wdth;
            
            for (int m = 0; m < 8; m++) // loop through 8-connected environment of the current edge point
            {
                int neighbourXPos = centreXPos + dX[m];
                int neighbourYPos = centreYPos + dY[m];
                if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds
                int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;
                int neighbourTag = cannyEdgeVector[neighbourIndex];
                if (neighbourTag == 3)
                {
                    cannyEdgeVector[neighbourIndex] = 6; // needed in case of cyclic path
                    connectedPoints.push_back(neighbourIndex);
                }
            }
        }
        
        // Run through all connected branches
        int numBranches = connectedPoints.size();
        for (int iBranch = 0; iBranch < numBranches; iBranch++)
        {
            int centreIndex = connectedPoints[iBranch];
            int vertexNeighbourIndex;
            
            // Create new branch, connect it with vertex and add first point
            branchProperties branchNew;
            branchNew.connectedVertices.resize(2); // each branch is connected to two vertices
            branchNew.connectedVertices[0] = iVertex;
            branchNew.pointIndices.push_back(centreIndex);
            
            // Scan through branch until new vertex is found
            bool FOUND_VERTEX = false;
            bool FOUND_POINT  = true;
            
            while (!FOUND_VERTEX && FOUND_POINT)
            {
                FOUND_POINT = false;
                vertexNeighbourIndex = -1;
                
                int centreXPos = centreIndex % mAOI.wdth;
                int centreYPos = (centreIndex - centreXPos) / mAOI.wdth;
                
                for (int m = 0; m < 8 && !FOUND_POINT; m++) // loop through 8-connected environment of the current edge point
                {
                    int neighbourXPos = centreXPos + dX[m];
                    int neighbourYPos = centreYPos + dY[m];
                    
                    if (neighbourXPos < 0 || neighbourXPos >= mAOI.wdth || neighbourYPos < 0 || neighbourYPos >= mAOI.hght) { continue; } // neighbour is out-of-bounds
                    
                    int neighbourIndex = mAOI.wdth * neighbourYPos + neighbourXPos;
                    
                    int neighbourTag = cannyEdgeVector[neighbourIndex];
                    
                    if (neighbourTag == 3) // branch continues
                    {
                        branchNew.pointIndices.push_back(neighbourIndex);
                        cannyEdgeVector[neighbourIndex] = 4; // new tag so that branch is not recorded again
                        centreIndex = neighbourIndex;
                        FOUND_POINT = true;
                    }
                    else if (neighbourTag == 5 || neighbourTag == 6) // check if vertex has been found
                    {
                        int jVertex = iVertex;
                        if (neighbourTag == 5) { jVertex = find2(vertexPointIndices, neighbourIndex); }
                        if (jVertex < numVertices)
                        {
                            int branchLength = branchNew.pointIndices.size();
                            if (jVertex != iVertex || branchLength > mDetectionParameters.windowLengthEdge)
                            {
                                branchNew.connectedVertices[1] = jVertex;
                                FOUND_VERTEX = true;
                                FOUND_POINT  = true;
                                continue;
                            }
                            vertexNeighbourIndex = jVertex;
                        }
                    }
                }
            }
            
            if (!FOUND_VERTEX) // no new vertex found
            {
                if (vertexNeighbourIndex >= 0) // cyclic path, add to existing vertex
                {
                    verticesAll[vertexNeighbourIndex].pointIndices.push_back(centreIndex);
                    vertexPointIndices[vertexNeighbourIndex].push_back(centreIndex);
                    branchNew.connectedVertices[1] = vertexNeighbourIndex;
                }
                else // create new vertex terminal
                {
                    vertexProperties vertexNew;
                    vertexNew.pointIndices.push_back(centreIndex);
                    vertexNew.tag = 1;
                    vertexPointIndices.push_back(vertexNew.pointIndices);
                    int jVertex = numVertices; // add to end
                    branchNew.connectedVertices[1] = jVertex;
                    verticesAll.push_back(vertexNew);
                    numVertices++;
                }
                
                branchNew.pointIndices.erase(branchNew.pointIndices.end() - 1);
                cannyEdgeVector[centreIndex] = 5;
            }
            
            branchesAll.push_back(branchNew);
        }
        
        for (int iEdgePoint = 0; iEdgePoint < numBranches; iEdgePoint++)
        { cannyEdgeVector[connectedPoints[iEdgePoint]] = 5; } // give normal tag now
    }
    
    // Record branches each vertex is connected to
    
    for (int iBranch = 0, numBranches = branchesAll.size(); iBranch < numBranches; iBranch++)
    {
        for (int iVertex = 0, vertexIndexOld = -1; iVertex < 2; iVertex++)
        {
            int vertexIndexNew = branchesAll[iBranch].connectedVertices[iVertex];
            if (vertexIndexOld != vertexIndexNew && vertexIndexNew < numVertices)
            { verticesAll[vertexIndexNew].connectedBranches.push_back(iBranch); }
            vertexIndexOld = vertexIndexNew;
        }
    }
    
    return branchesAll;
}

std::vector<std::vector<int>> depthFirstSearch(std::vector<vertexProperties>& vVertexPropertiesAll, std::vector<branchProperties>& vBranchPropertiesAll, std::vector<int>& pathBranchesRoot, std::vector<int>& verticesCheckedAll, std::vector<int>& verticesChecked, std::vector<int>& branchesChecked, int vertexIndex)
{
    // Depth-first search
    std::vector<std::vector<int>> pathBranchesAll;
    
    vertexProperties mVertex = vVertexPropertiesAll[vertexIndex];
    int numBranches = mVertex.connectedBranches.size();
    
    for (int iBranch = 0; iBranch < numBranches; iBranch++)
    {
        int branchIndexNew = mVertex.connectedBranches[iBranch];
        
        if (branchesChecked[branchIndexNew] == 0)
        {
            std::vector<int> pathBranchesRootNew = pathBranchesRoot;
            pathBranchesRootNew.push_back(branchIndexNew);
            
            std::vector<int> branchesCheckedNew = branchesChecked;
            branchesCheckedNew[branchIndexNew]  = 1;
            
            std::vector<int> verticesCheckedNew = verticesChecked;
            verticesCheckedNew[vertexIndex]     = 1;
            
            // Grab other vertex edge is attached to
            
            std::vector<int> connectedVertices = vBranchPropertiesAll[branchIndexNew].connectedVertices;
            
            int vertexIndexNew;
            if (connectedVertices[0] == vertexIndex) { vertexIndexNew = connectedVertices[1]; }
            else                                     { vertexIndexNew = connectedVertices[0]; }
            
            // Find all paths from next vertex
            
            if (verticesCheckedAll[vertexIndexNew] == 0) // don't add branch
            {
                if (verticesCheckedNew[vertexIndexNew] == 0) // add branch but stop at vertex
                {
                    std::vector<std::vector<int>> pathBranchesNew = depthFirstSearch(vVertexPropertiesAll, vBranchPropertiesAll, pathBranchesRootNew, verticesCheckedAll, verticesCheckedNew, branchesCheckedNew, vertexIndexNew);
                    if (pathBranchesNew.size() > 0) { pathBranchesAll.insert(std::end(pathBranchesAll), std::begin(pathBranchesNew), std::end(pathBranchesNew)); }
                }
                
                pathBranchesAll.push_back(pathBranchesRootNew); // always add root path
            }
        }
    }
    
    return pathBranchesAll;
}

bool findPreferredPath(int& acceptedPathIndex, const std::vector<branchProperties>& vBranchProperties, const std::vector<std::vector<int>>& pathsAll, const double& prediction, const double& lowerLimit, const double& upperLimit, const bool& CYCLIC)
{
    int numPaths = pathsAll.size();
    std::vector<int> pathLengths(numPaths);
    
    for (int iPath = 0; iPath < numPaths; iPath++)
    {
        int lengthTotal = 0;
        std::vector<int> pathCurrent = pathsAll[iPath];
        int numPathBranches = pathCurrent.size();
        
        for (int iBranch = 0; iBranch < numPathBranches; iBranch++)
        {
            lengthTotal += vBranchProperties[pathCurrent[iBranch]].length;
        }
        
        if (!CYCLIC || (lengthTotal > lowerLimit && lengthTotal < upperLimit))
        {
            pathLengths[iPath] = lengthTotal;
        }
    }
    
    int numPathsNew = pathLengths.size();
    
    if (pathLengths.size() > 0)
    {
        acceptedPathIndex = 0;
        
        if (numPathsNew > 1) // If multiple paths, choose closest to prediction
        {
            std::vector<double> lengthError(numPathsNew);
            
            for (int iPath = 0; iPath < numPathsNew; iPath++)
            {
                lengthError[iPath] = std::abs(pathLengths[iPath] - prediction);
            }
            
            acceptedPathIndex = std::distance(lengthError.begin(), std::min_element(lengthError.begin(), lengthError.end()));
        }
        
        return true;
    }
    else { return false; }
}

std::vector<int> processGraphTree(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, std::vector<vertexProperties>& vVertexPropertiesAll, std::vector<branchProperties>& vBranchPropertiesAll, const AOIProperties& mAOI)
{
    int numBranchesAll = vBranchPropertiesAll.size();
    int numVerticesAll = vVertexPropertiesAll.size();
    
    // Perform depth-first search for each vertex
    
    std::vector<std::vector<int>> pathsAll;
    std::vector<int> cyclicPaths;
    
    std::vector<int> verticesCheckedAll(numVerticesAll, 0);
    
    // Do terminal vertices first
    
    std::vector<int> terminalVertices;
    for (int iVertex = 0; iVertex < numVerticesAll; iVertex++)
    { if (vVertexPropertiesAll[iVertex].tag == 1) { terminalVertices.push_back(iVertex); }}
    int numTerminals = terminalVertices.size();
    
    for (int iVertex = 0; iVertex < numTerminals; iVertex++) // loop through all (starting) vertices
    {
        int jVertex = terminalVertices[iVertex];
        if (vVertexPropertiesAll[jVertex].connectedBranches.size() > 0) // ignore isolated nodes
        {
            std::vector<int> branchesChecked(numBranchesAll, 0); // no branches checked
            std::vector<int> verticesChecked(numVerticesAll, 0); // no branches checked
            std::vector<int> pathBranchesRoot; // start with no branches
            std::vector<std::vector<int>> pathBranchesAll = depthFirstSearch(vVertexPropertiesAll, vBranchPropertiesAll, pathBranchesRoot, verticesCheckedAll, verticesChecked, branchesChecked, jVertex);
            pathsAll.insert(std::end(pathsAll), std::begin(pathBranchesAll), std::end(pathBranchesAll));
        }
        
        verticesCheckedAll[jVertex] = 1; // never end with a vertex terminal that we already started with (redundant)
    }
    
    // Then do branch vertices
    
    for (int iVertex = 0; iVertex < numVerticesAll; iVertex++) // loop through all (starting) vertices
    {
        if (verticesCheckedAll[iVertex] == 0)
        {
            if (vVertexPropertiesAll[iVertex].connectedBranches.size() > 0) // ignore isolated nodes
            {
                std::vector<int> branchesChecked(numBranchesAll, 0); // no branches checked
                std::vector<int> verticesChecked(numVerticesAll, 0); // no branches checked
                std::vector<int> pathBranchesRoot; // start with no branches
                std::vector<std::vector<int>> pathBranchesAll = depthFirstSearch(vVertexPropertiesAll, vBranchPropertiesAll, pathBranchesRoot, verticesCheckedAll, verticesChecked, branchesChecked, iVertex);
                pathsAll.insert(std::end(pathsAll), std::begin(pathBranchesAll), std::end(pathBranchesAll));

                // Check for cyclic path
                if (vVertexPropertiesAll[iVertex].tag == 2) // only internal vertices can have cyclic path
                {
                    int numPaths = pathBranchesAll.size();
                    if (numPaths > 1)
                    {
                        for (int iPath = 0; iPath < numPaths; iPath++)
                        {
                            std::vector<int> connectedBranches = vVertexPropertiesAll[iVertex].connectedBranches; // all branches vertex connects
                            std::vector<int>::iterator itr = find(connectedBranches.begin(), connectedBranches.end(), pathBranchesAll[iPath].back()); // check if last added branch connects to starting vertex
                            if (itr != connectedBranches.end()) { cyclicPaths.push_back(iPath); }
                        }
                    }
                }
            }
        }
    }
    
    // Find optimal path
    
    double circumferencePrediction = mDetectionVariables.predictedCircumference;
    
    // Get length of each branch
    
    for (int iBranch = 0; iBranch < numBranchesAll; iBranch++)
    {
        vBranchPropertiesAll[iBranch].length = vBranchPropertiesAll[iBranch].pointIndices.size();
    }
    
    std::vector<int> pathBranchIndices;
    
    // First consider cyclic paths
    
    int numCyclicPaths = cyclicPaths.size();
    
    if (numCyclicPaths > 0)
    {
        double circumferenceLowerLimit = circumferencePrediction * (1 - mDetectionVariables.thresholdChangeCircumference);
        double circumferenceUpperLimit = circumferencePrediction * (1 + mDetectionVariables.thresholdChangeCircumference);
        
        if (circumferenceLowerLimit < mDetectionParameters.circumferenceMin) { circumferenceLowerLimit = mDetectionParameters.circumferenceMin; }
        if (circumferenceUpperLimit > mDetectionParameters.circumferenceMax) { circumferenceUpperLimit = mDetectionParameters.circumferenceMax; }
        
        std::vector<std::vector<int>> cyclicPathsAll(numCyclicPaths);
        for (int iPath = 0; iPath < numCyclicPaths; iPath++) { cyclicPathsAll[iPath] = pathsAll[cyclicPaths[iPath]]; }
        
        int acceptedPathIndex;
        if (findPreferredPath(acceptedPathIndex, vBranchPropertiesAll, cyclicPathsAll, circumferencePrediction, circumferenceLowerLimit, circumferenceUpperLimit, true))
        { pathBranchIndices = cyclicPathsAll[acceptedPathIndex]; }
    }
    
    int numBranchesPath = pathBranchIndices.size();
    
    // If no cyclic path found, then take path closest to prediction length
    
    if (numBranchesPath == 0)
    {
        int acceptedPathIndex;
        findPreferredPath(acceptedPathIndex, vBranchPropertiesAll, pathsAll, circumferencePrediction);
        pathBranchIndices = pathsAll[acceptedPathIndex];
        numBranchesPath   = pathBranchIndices.size();
    }
    
    // Grab branches for path
    
    std::vector<branchProperties> vBranchPropertiesPath(numBranchesPath);
    
    for (int iBranch = 0; iBranch < numBranchesPath; iBranch++)
    {
        vBranchPropertiesPath[iBranch] = vBranchPropertiesAll[pathBranchIndices[iBranch]];
    }
    
    // Add vertices to path, and
    // reverse branches that are not properly aligned
    
    if (numBranchesPath > 1)
    {
        for (int iBranch = 0; iBranch < numBranchesPath - 1; iBranch++)
        {
            branchProperties mBranchProperties_1 = vBranchPropertiesPath[iBranch];
            branchProperties mBranchProperties_2 = vBranchPropertiesPath[iBranch + 1];
            
            for (int iVertex = 0; iVertex < 2; iVertex++)
            {
                int vertex_1 = mBranchProperties_1.connectedVertices[iVertex];
                
                for (int jVertex = 0; jVertex < 2; jVertex++)
                {
                    int vertex_2 = mBranchProperties_2.connectedVertices[jVertex];
                    if (vertex_1 == vertex_2) // last vertex of current branch should be equal to first vertex of next branch
                    {
                        if (iVertex == 0) { std::reverse(vBranchPropertiesPath[iBranch].pointIndices.begin(), vBranchPropertiesPath[iBranch].pointIndices.end()); }
                        vBranchPropertiesPath[iBranch].pointIndices.insert(vBranchPropertiesPath[iBranch].pointIndices.begin(),
                                                                           vVertexPropertiesAll[vertex_1].pointIndices.begin(),
                                                                           vVertexPropertiesAll[vertex_1].pointIndices.end());
                        
                        if (iBranch == 0)
                        {
                            int vertexIndex = mBranchProperties_1.connectedVertices[(iVertex + 1) % 2];
                            vBranchPropertiesPath[iBranch].pointIndices.insert(vBranchPropertiesPath[iBranch].pointIndices.begin(),
                                                                               vVertexPropertiesAll[vertexIndex].pointIndices.begin(),
                                                                               vVertexPropertiesAll[vertexIndex].pointIndices.end());
                        }
                        
                        if (iBranch == numBranchesPath - 2)
                        {
                            int vertexIndex = mBranchProperties_2.connectedVertices[(jVertex + 1) % 2];
                            vBranchPropertiesPath[iBranch + 1].pointIndices.insert(vBranchPropertiesPath[iBranch + 1].pointIndices.end(),
                                    vVertexPropertiesAll[vertexIndex].pointIndices.begin(),
                                    vVertexPropertiesAll[vertexIndex].pointIndices.end());
                        }
                        
                        break;
                    }
                }
            }
        }
    }
    else
    {
        std::vector<int> dZ = {-1, -mAOI.wdth - 1, -mAOI.wdth, -mAOI.wdth + 1, 1, mAOI.wdth + 1,mAOI.wdth, mAOI.wdth - 1};
        
        branchProperties mBranchPropertiesPath = vBranchPropertiesPath[0];
        int numBranchPoints = mBranchPropertiesPath.pointIndices.size();
        
        int vertex_1 = mBranchPropertiesPath.connectedVertices[0];
        int vertex_2 = mBranchPropertiesPath.connectedVertices[1];
        std::vector<int> vertexPoints_1  = vVertexPropertiesAll[vertex_1].pointIndices;
        std::vector<int> vertexPoints_2  = vVertexPropertiesAll[vertex_2].pointIndices;
        
        if (numBranchPoints > 0)
        {
            std::vector<int> vertexPointsBegin, vertexPointsEnd;
            int pointIndex = mBranchPropertiesPath.pointIndices[0];
            for (int m = 0; m < 8; m++)
            {
                int neighbourIndex = pointIndex + dZ[m];
                
                std::vector<int>::iterator itr_1 = find(vertexPoints_1.begin(), vertexPoints_1.end(), neighbourIndex); // check if index has already been stored
                if (itr_1 == vertexPoints_1.end())
                {
                    vertexPointsBegin = vertexPoints_1;
                    vertexPointsEnd   = vertexPoints_2;
                    break;
                }
                
                std::vector<int>::iterator itr_2 = find(vertexPoints_2.begin(), vertexPoints_2.end(), neighbourIndex); // check if index has already been stored
                if (itr_2 == vertexPoints_2.end())
                {
                    vertexPointsBegin = vertexPoints_2;
                    vertexPointsEnd   = vertexPoints_1;
                    break;
                }
            }
            
            mBranchPropertiesPath.pointIndices.insert(mBranchPropertiesPath.pointIndices.begin(), vertexPointsBegin.begin(), vertexPointsBegin.end());
            mBranchPropertiesPath.pointIndices.insert(mBranchPropertiesPath.pointIndices.end(),   vertexPointsEnd.begin(),   vertexPointsEnd.end());
        }
        else
        {
            mBranchPropertiesPath.pointIndices.insert(mBranchPropertiesPath.pointIndices.begin(), vertexPoints_1.begin(), vertexPoints_1.end());
            mBranchPropertiesPath.pointIndices.insert(mBranchPropertiesPath.pointIndices.end(),   vertexPoints_2.begin(), vertexPoints_2.end());
        }
        
        vBranchPropertiesPath[0] = mBranchPropertiesPath; // update
    }
    
    // All path points
    
    std::vector<int> pathPoints;
    
    for (int iBranch = 0; iBranch < numBranchesPath; iBranch++)
    {
        branchProperties mBranchProperties = vBranchPropertiesPath[iBranch];
        pathPoints.insert(std::end(pathPoints), std::begin(mBranchProperties.pointIndices), std::end(mBranchProperties.pointIndices));
    }
    
    return pathPoints;
}

std::vector<edgeProperties> edgeSelection(const detectionVariables& mDetectionVariables, detectionParameters& mDetectionParameters, std::vector<int>& cannyEdgeVector, AOIProperties mAOI)
{
    std::vector<edgeProperties> vEdgePropertiesAll; // new structure containing length and indices of all selected edges
    
    std::vector<int> startIndicesRaw = findEdges(mDetectionVariables, mDetectionParameters, cannyEdgeVector, mAOI);
    int numOrigins = startIndicesRaw.size();

    int numEdges = 0;

    do
    {
        std::vector<std::vector<int>> vAllIndices;
        for (int iEdge = 0; iEdge < numOrigins; iEdge++)
        {
            int startIndex = startIndicesRaw[iEdge];
            if (cannyEdgeVector[startIndex] == 1)
            {
                std::vector<int> edgePointIndices = findEdgePoints(mDetectionParameters, cannyEdgeVector, mAOI, startIndex); // tag all edge points of found edges
                vAllIndices.push_back(edgePointIndices);
            }
        }
        
        numEdges = vAllIndices.size();
        
        for (int iEdge = 0; iEdge < numEdges; iEdge++)
        {
            std::vector<int> allIndices = vAllIndices[iEdge];
            int numEdgePointsTotal = allIndices.size();
            if (numEdgePointsTotal > mDetectionParameters.windowLengthEdge)
            {
                int startIndex = allIndices.back();
                
                // Find all vertices and all connected branches (i.e. obtain graph tree)
                
                std::vector<vertexProperties> vVertexProperties = findGraphVertices(cannyEdgeVector, mAOI, startIndex);
                std::vector<branchProperties> vBranchProperties = findGraphBranches(mDetectionParameters, vVertexProperties, cannyEdgeVector, mAOI);
                
                // Find preferred path:
                // Cyclic path that resembles pupil outline the most,
                // otherwise take path closest to circumference prediction
                
                int numBranches = vBranchProperties.size();
                int numVertices = vVertexProperties.size();
                
                if (numBranches > 0 && numVertices > 0)
                {
                    std::vector<int> pathIndices = processGraphTree(mDetectionVariables, mDetectionParameters, vVertexProperties, vBranchProperties, mAOI);
                    
                    // Give points in optimal path a new tag
                    
                    for (int iEdgePoint = 0, edgeSize = pathIndices.size(); iEdgePoint < edgeSize; iEdgePoint++)
                    { cannyEdgeVector[pathIndices[iEdgePoint]] = 7; }
                    
                    edgeProperties mEdgeProperties;
                    mEdgeProperties.pointIndices = pathIndices;
                    
                    vEdgePropertiesAll.push_back(mEdgeProperties);
                    
                    // Remove tag from points that have been tagged before, but not included in final path
                    
                    for (int iEdgePoint = 0; iEdgePoint < numEdgePointsTotal; iEdgePoint++)
                    {
                        int edgePointIndex  = allIndices[iEdgePoint];
                        int edgePointTag    = cannyEdgeVector[edgePointIndex];
                        if (edgePointTag >= 2 && edgePointTag <= 6)
                        { cannyEdgeVector[edgePointIndex] = 1; }
                    }
                }
            }
        }
    } while (numEdges > 1);
    
    return vEdgePropertiesAll;
}

void calculateCurvatureLimits(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, double& curvatureUpperLimit, double& curvatureLowerLimit)
{
    // Calculate curvature limits
    
    std::vector<double> circumferences(2);
    std::vector<double>   aspectRatios(2);
    
    circumferences[0] = mDetectionVariables.predictedCircumference + (mDetectionVariables.predictedCircumference * mDetectionVariables.thresholdChangeCircumference);
    circumferences[1] = mDetectionVariables.predictedCircumference - (mDetectionVariables.predictedCircumference * mDetectionVariables.thresholdChangeCircumference);

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
    
    curvatureUpperLimit = *std::max_element(std::begin(curvaturesMax), std::end(curvaturesMax)) + (M_PI * mDetectionParameters.curvatureOffset / 180);
    curvatureLowerLimit = *std::min_element(std::begin(curvaturesMin), std::end(curvaturesMin)) - (M_PI * mDetectionParameters.curvatureOffset / 180);
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
        
        curvatures[iEdgePoint] = vectorAngle / mDetectionParameters.windowLengthEdge;
        
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
    else // ignore curvature
    {
        curvatureAvg = std::numeric_limits<double>::quiet_NaN();
        curvatureMax = std::numeric_limits<double>::quiet_NaN();
        curvatureMin = std::numeric_limits<double>::quiet_NaN();
    }
    
    mEdgeProperties.curvature    = curvatureAvg;
    mEdgeProperties.curvatureMax = curvatureMax;
    mEdgeProperties.curvatureMin = curvatureMin;
}

std::vector<edgeProperties> edgeSegmentationLength(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const edgeProperties& mEdgeProperties)
{
    // This functions cuts edge terminals to make the edge shorter if the edge is significantly longer than predicted
    
    int edgeSize = mEdgeProperties.curvatures.size();
    
    // find breakpoints based on length thresholding
    
    std::vector<int> breakPoints; // position of breakpoints
    breakPoints.push_back(-1); // add first point (+ 1 is added later)
    
    double lengthDifference = mEdgeProperties.length - mDetectionVariables.predictedCircumference;
    
    if (lengthDifference > mDetectionParameters.windowLengthEdge && edgeSize > 2 * lengthDifference + mDetectionParameters.windowLengthEdge) // should be enough space between breakpoints
    {
        breakPoints.push_back(lengthDifference);
        breakPoints.push_back(edgeSize - lengthDifference - 1);
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
            mEdgePropertiesNew.radiusVar = calculateVariance(mEdgeProperties.radii) / mDetectionVariables.predictedCircumference;
            mEdgePropertiesNew.length    = mEdgePropertiesNew.pointIndices.size();
            
            calculateCurvatureStats(mDetectionParameters, mEdgePropertiesNew);
            
            vEdgeProperties[iBreakPoint] = mEdgePropertiesNew;
        }
        
        // only remove one of the two edge terminals - re-attach other one
        // compare edge characteristics of terminals with middle section
        // 0 = start terminal
        // 1 = main section
        // 2 = end terminal
        
        std::vector<double> featureValues_1(6); // for start terminal
        featureValues_1[0] = std::abs(vEdgeProperties[0].radius     - vEdgeProperties[1].radius)                  / std::max(vEdgeProperties[0].radius, vEdgeProperties[1].radius);
        featureValues_1[1] = std::abs(vEdgeProperties[0].length     - mDetectionVariables.predictedCircumference) / std::max(vEdgeProperties[0].length, mDetectionVariables.predictedCircumference);
        featureValues_1[2] = std::abs(vEdgeProperties[0].curvature  - vEdgeProperties[1].curvature);
        featureValues_1[3] = std::abs(vEdgeProperties[0].intensity  - vEdgeProperties[1].intensity);
        featureValues_1[4] = std::abs(vEdgeProperties[0].gradient   - vEdgeProperties[1].gradient);
        featureValues_1[5] = std::abs(vEdgeProperties[0].radiusVar  - vEdgeProperties[1].radiusVar);
        
        std::vector<double> featureValues_2(6); // for end terminal
        featureValues_2[0] = std::abs(vEdgeProperties[2].radius     - vEdgeProperties[1].radius)                  / std::max(vEdgeProperties[2].radius, vEdgeProperties[1].radius);
        featureValues_2[1] = std::abs(vEdgeProperties[2].length     - mDetectionVariables.predictedCircumference) / std::max(vEdgeProperties[2].length, mDetectionVariables.predictedCircumference);
        featureValues_2[2] = std::abs(vEdgeProperties[2].curvature  - vEdgeProperties[1].curvature);
        featureValues_2[3] = std::abs(vEdgeProperties[2].intensity  - vEdgeProperties[1].intensity);
        featureValues_2[4] = std::abs(vEdgeProperties[2].gradient   - vEdgeProperties[1].gradient);
        featureValues_2[5] = std::abs(vEdgeProperties[2].radiusVar  - vEdgeProperties[1].radiusVar);
        
        double scoreTotal_1 = calculateScoreTotal(mDetectionVariables, featureValues_1, false, false);
        double scoreTotal_2 = calculateScoreTotal(mDetectionVariables, featureValues_2, false, false);
        
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
        
        if (offsetXPos < 0 || offsetXPos >= mAOI.wdth || offsetYPos < 0 || offsetYPos >= mAOI.hght)
        {       edgeIntensities[iEdgePoint] = (int) ptr_img[edgePointXPos + edgePointYPos * mAOI.wdth]; }
        else {  edgeIntensities[iEdgePoint] = (int) ptr_img[   offsetXPos +    offsetYPos * mAOI.wdth]; }
    }
    
    return edgeIntensities;
}

std::vector<edgeProperties> edgeSegmentationScore(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const edgeProperties& mEdgeProperties)
{
    // start from end and move through edge with window
    // if window score average is significantly below central average --> remove
    
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
            std::vector<double> scoreDifferenceVector;
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
                    mEdgePropertiesNew.radiusVar = calculateVariance  (mEdgePropertiesNew.radii) / mDetectionVariables.predictedCircumference;
                    mEdgePropertiesNew.intensity = calculateMeanInt   (mEdgePropertiesNew.intensities);
                    mEdgePropertiesNew.gradient  = calculateMeanInt   (mEdgePropertiesNew.gradients);
                    mEdgePropertiesNew.length    = mEdgePropertiesNew.pointIndices.size();
                    calculateCurvatureStats(mDetectionParameters, mEdgePropertiesNew); // get average curvature
                    
                    vEdgePropertiesTemp[i] = mEdgePropertiesNew;
                }
                
                // Calculate score difference
                
                double predictedRadius = mDetectionVariables.predictedCircumference / (2 * M_PI);
                
                std::vector<double> featureValues_1(6);
                featureValues_1[0] = std::abs(vEdgePropertiesTemp[0].radius    - predictedRadius)                            / std::max(vEdgePropertiesTemp[0].radius, predictedRadius);
                featureValues_1[1] = std::abs(vEdgePropertiesTemp[0].length    - mDetectionVariables.predictedCircumference) / std::max(vEdgePropertiesTemp[0].length, mDetectionVariables.predictedCircumference);
                featureValues_1[2] = std::abs(vEdgePropertiesTemp[0].curvature - mDetectionVariables.predictedCurvature);
                featureValues_1[3] = std::abs(vEdgePropertiesTemp[0].intensity - mDetectionVariables.predictedIntensity);
                featureValues_1[4] = std::abs(vEdgePropertiesTemp[0].gradient  - mDetectionVariables.predictedGradient);
                featureValues_1[5] = std::abs(vEdgePropertiesTemp[0].radiusVar);
                
                std::vector<double> featureValues_2(6);
                featureValues_2[0] = std::abs(vEdgePropertiesTemp[1].radius    - predictedRadius)                            / std::max(vEdgePropertiesTemp[1].radius, predictedRadius);
                featureValues_2[1] = std::abs(vEdgePropertiesTemp[1].length    - mDetectionVariables.predictedCircumference) / std::max(vEdgePropertiesTemp[1].length, mDetectionVariables.predictedCircumference);
                featureValues_2[2] = std::abs(vEdgePropertiesTemp[1].curvature - mDetectionVariables.predictedCurvature);
                featureValues_2[3] = std::abs(vEdgePropertiesTemp[1].intensity - mDetectionVariables.predictedIntensity);
                featureValues_2[4] = std::abs(vEdgePropertiesTemp[1].gradient  - mDetectionVariables.predictedGradient);
                featureValues_2[5] = std::abs(vEdgePropertiesTemp[1].radiusVar);
                
                double score_1 = calculateScoreTotal(mDetectionVariables, featureValues_1, false, true);
                double score_2 = calculateScoreTotal(mDetectionVariables, featureValues_2, false, true);
                
                if (score_1 < mDetectionVariables.thresholdScore && score_2 < mDetectionVariables.thresholdScore) { break; } // ignore probable non-pupil edges
                
                double scoreDifference = std::abs(score_1 - score_2);
                
                if (scoreDifference < mDetectionParameters.thresholdScoreDiffEdge)
                {
                    double nSteps = (mDetectionParameters.thresholdScoreDiffEdge - scoreDifference) / (scoreStepSize * mDetectionParameters.thresholdScoreDiffEdge);
                    iEdgePoint    = iEdgePoint + nSteps * mDetectionParameters.windowLengthEdge;// skip next division(s)
                }
                else
                {
                    scoreDifferenceVector.push_back(std::abs(scoreDifference));
                    breakPointsAll.push_back(iEdgePoint);
                    iBreakPoint++;
                }
            } while (!BREAK_LOOP);
            
            // Find point of maximum score difference
            
            auto itr = std::max_element(scoreDifferenceVector.begin(), scoreDifferenceVector.end());
            double scoreDifferenceMax = *itr;
            int indexMax = std::distance(scoreDifferenceVector.begin(), itr);
            
            std::vector<edgeProperties> vEdgePropertiesTemp;
            
            if (scoreDifferenceMax > mDetectionParameters.thresholdScoreDiffEdge)
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
                cannyEdgeVector[neighbourIndex] = 7; // ... tag it and ...
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
    
    if (mDetectionVariables.certaintyFeatures > certaintyThreshold || mDetectionVariables.certaintyPosition > certaintyThreshold) // need to have enough certainty
    {
        for (int iEdge = 0; iEdge < numEdges; iEdge++)
        {
            std::vector<double> featureValues(6);
            
            double radius = vEdgePropertiesAll[iEdge].radius;
            double length = vEdgePropertiesAll[iEdge].length;
            
            double radiusPredicted = mDetectionVariables.predictedCircumference / (2 * M_PI);
            
            featureValues[0] = std::abs(radius - radiusPredicted) / std::max(radius, radiusPredicted);
            featureValues[1] = std::abs(length - mDetectionVariables.predictedCircumference) / std::max(length, mDetectionVariables.predictedCircumference);
            featureValues[2] = std::abs(vEdgePropertiesAll[iEdge].curvature - mDetectionVariables.predictedCurvature);
            featureValues[3] = std::abs(vEdgePropertiesAll[iEdge].intensity - mDetectionVariables.predictedIntensity);
            featureValues[4] = std::abs(vEdgePropertiesAll[iEdge].gradient  - mDetectionVariables.predictedGradient);
            featureValues[5] = vEdgePropertiesAll[iEdge].radiusVar;
            
            totalScores[iEdge] = calculateScoreTotal(mDetectionVariables, featureValues, true, true);
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

std::vector<edgeProperties> edgeCollectionFilter(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const std::vector<edgeProperties>& vEdgePropertiesAll, const AOIProperties& mAOI)
{
    std::vector<edgeProperties> vEdgeProperties; // properties of edge collections
    
    int numEdgesTotal = vEdgePropertiesAll.size(); // total number of edges
    
    // First collect all ellipse fits
    
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
            
            std::vector<int> edgePointIndices; // vector containing all point indices for fit
            edgePointIndices.reserve(edgeSetSize); // preallocate memory
            
            for (int iEdge = 0; iEdge < combiNumEdges; iEdge++)
            { edgePointIndices.insert(edgePointIndices.end(), combiEdgePoints[iEdge].begin(), combiEdgePoints[iEdge].end()); }
            
            // Calculate range
            
            std::vector<int> edgeXPositions(edgeSetSize);
            std::vector<int> edgeYPositions(edgeSetSize);
            
            for (int iEdgePoint = 0; iEdgePoint < edgeSetSize; iEdgePoint++)
            {
                int edgePointIndex = edgePointIndices[iEdgePoint];
                edgeXPositions[iEdgePoint] = edgePointIndex % mAOI.wdth;
                edgeYPositions[iEdgePoint] = (edgePointIndex - edgeXPositions[iEdgePoint]) / mAOI.wdth;
            }
            
            int XPosMin = *std::min_element(std::begin(edgeXPositions), std::end(edgeXPositions));
            int XPosMax = *std::max_element(std::begin(edgeXPositions), std::end(edgeXPositions));
            int YPosMin = *std::min_element(std::begin(edgeYPositions), std::end(edgeYPositions));
            int YPosMax = *std::max_element(std::begin(edgeYPositions), std::end(edgeYPositions));
            
            double combiWdth = 0.5 * (XPosMax - XPosMin);
            double combiHght = 0.5 * (YPosMax - YPosMin);

            if (combiWdth < fitMinRange * mDetectionVariables.predictedWidth || combiHght < fitMinRange * mDetectionVariables.predictedHeight) { continue; }
            
            // Record edge collection
            
            edgeProperties mEdgeProperties;
            mEdgeProperties.edgeIndices  = combiEdgeIndices;
            mEdgeProperties.pointIndices = edgePointIndices;
            mEdgeProperties.length       = edgeSetLength;
            mEdgeProperties.size         = edgeSetSize;
            vEdgeProperties.push_back(mEdgeProperties);
            
        } while (std::next_permutation(edgeCombination.begin(), edgeCombination.end()));
    }
    
    // No more than fixed number of fits
    
    int numFits = vEdgeProperties.size();
    if (numFits > mDetectionParameters.fitMaximum)
    {
        std::vector<edgeProperties> vEdgePropertiesTemp(mDetectionParameters.fitMaximum);
        std::vector<double> lengthErrorUnsorted(numFits);
        
        for (int iFit = 0; iFit < numFits; iFit++)
        { lengthErrorUnsorted[iFit] = std::abs(vEdgeProperties[iFit].length - mDetectionVariables.predictedCircumference); }
        
        std::vector<double> lengthErrorSorted = lengthErrorUnsorted;
        std::sort(lengthErrorSorted.begin(), lengthErrorSorted.end());
        
        for (int iFit = 0; iFit < mDetectionParameters.fitMaximum; iFit++)
        {
            for (int jFit = 0; jFit < numFits; jFit++)
            {
                if (lengthErrorSorted[iFit] == lengthErrorUnsorted[jFit])
                {
                    vEdgePropertiesTemp[iFit] = vEdgePropertiesTemp[jFit];
                    lengthErrorUnsorted[jFit] = std::numeric_limits<double>::max();
                    break;
                }
            }
        }
        
        vEdgeProperties = vEdgePropertiesTemp;
    }
    
    return vEdgeProperties;
}

std::vector<ellipseProperties> ellipseFitting(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, const std::vector<edgeProperties>& vEdgePropertiesAll, AOIProperties mAOI)
{
    std::vector<ellipseProperties> vEllipsePropertiesAll; // vector to record information for each accepted ellipse fit
    
    int numEdgesTotal = vEdgePropertiesAll.size();
    
    for (int iEdge = 0; iEdge < numEdgesTotal; iEdge++)
    {
        std::vector<int> edgeIndices = vEdgePropertiesAll[iEdge].edgeIndices;
        double edgeSetLength = vEdgePropertiesAll[iEdge].length;
        double edgeSetSize   = vEdgePropertiesAll[iEdge].size;

        ellipseProperties mEllipseProperties = fitEllipse(edgeIndices, edgeSetSize, mAOI);
        
        if (!mEllipseProperties.DETECTED) { continue; } // error
        
        if (edgeSetLength < fitMinEdgeLength * mEllipseProperties.circumference) { continue; } // minimum number of edge points required
        
        // Absolute size and shape filter
        
        if (mEllipseProperties.circumference < mDetectionParameters.circumferenceMin) { continue; } // no small ellipse
        if (mEllipseProperties.aspectRatio   < mDetectionParameters.aspectRatioMin) { continue; } // no extreme deviations from circular shape
        
        // Size and shape combination filter
        
        double circumferenceUpperLimit = aspectRatioSlope * (mEllipseProperties.aspectRatio - 1) + mDetectionParameters.circumferenceMax;
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
        mEllipseProperties.edgeIndices = edgeIndices;
        mEllipseProperties.edgeLength  = edgeSetLength;
        vEllipsePropertiesAll.push_back(mEllipseProperties);
    }
    
    return vEllipsePropertiesAll;
}

std::vector<int> ellipseFitFilter(const detectionVariables& mDetectionVariables, const detectionParameters& mDetectionParameters, std::vector<ellipseProperties> vEllipseProperties)
{
    static const double weightAspectRatio   = 0.815;
    static const double weightCircumference = 0.415;
    static const double weightLength        = 1.325;
    static const double weightError         = 0.239;
    static const double weightAngle         = 1.014;
    static const double weightAngleFactor   = 0.787;

    static const double parametersAspectRatio   = 0.70922;
    static const double parametersCircumference = 0.74258;
    static const double parametersLength        = 0.32519;
    static const double parametersError         = 0.28342;
    static const double parametersAngle         = 0.31997;
    
    int numFits = vEllipseProperties.size();
    
    std::vector<double> scoreFits(numFits);
    
    for (int iFit = 0; iFit < numFits; iFit++)
    {
        ellipseProperties mEllipseProperties = vEllipseProperties[iFit];
        
        double errorAspectRatio   = std::abs(mEllipseProperties.aspectRatio - mDetectionVariables.predictedAspectRatio);
        double errorAngle         = std::abs(mEllipseProperties.angle       - mDetectionVariables.predictedAngle);
        double ratioCircumference = mEllipseProperties.circumference / mDetectionVariables.predictedCircumference;
        double ratioLength        = mEllipseProperties.edgeLength    / mDetectionVariables.predictedCircumference;
        double fitError           = mEllipseProperties.fitError;
        
        double factorAngleFunction = 1 - weightAngleFactor * mEllipseProperties.aspectRatio;
        
        double factorAngle         = mDetectionVariables.certaintyFeatures * weightAngle         * factorAngleFunction;
        double factorAspectRatio   = mDetectionVariables.certaintyFeatures * weightAspectRatio;
        double factorCircumference = mDetectionVariables.certaintyFeatures * weightCircumference;
        double factorLength        = mDetectionVariables.certaintyFeatures * weightLength;
        double factorFitError      =                           weightError;
        
        // Calculate scores
        
        double scoreAngle         = factorAngle         * gaussian(errorAngle,         parametersAngle         );
        double scoreAspectRatio   = factorAspectRatio   * gaussian(errorAspectRatio,   parametersAspectRatio   );
        double scoreCircumference = factorCircumference * gaussian(ratioCircumference, parametersCircumference );
        double scoreLength        = factorLength        * gaussian(ratioLength,        parametersLength        );
        double scoreFitError      = factorFitError      * gaussian(fitError,           parametersError         );
        
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
    
    // Update position prediction using Haar feature detector response




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
    
    std::vector<int> cannyEdgesOriginal  = cannyConversion(imageCannyEdges, cannyAOI); // convert to binary vector
    std::vector<int> cannyEdgesSharpened = cannyEdgesOriginal;
    std::vector<int> edgePointsOriginal  = getEdgeIndices(cannyEdgesOriginal, 1);
    std::vector<int> edgePointsSharpened;
    edgePointsSharpened = sharpenEdges_1(cannyEdgesSharpened, edgePointsOriginal,  mAOI);
    edgePointsSharpened = sharpenEdges_2(cannyEdgesSharpened, edgePointsSharpened, mAOI);
    
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
    
    if (mDetectionVariables.certaintyPosition > certaintyThreshold || mDetectionVariables.certaintyFeatures > certaintyThreshold)
    {
        std::vector<edgeProperties> vEdgePropertiesNew;
        
        for (int iEdge = 0, numEdges = vEdgePropertiesAll.size(); iEdge < numEdges; iEdge++)
        {
            edgeProperties mEdgeProperties    = vEdgePropertiesAll[iEdge];
            
            std::vector<edgeProperties> vEdgePropertiesTemp = edgeSegmentationScore(mDetectionVariables, mDetectionParameters, mEdgeProperties);
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
            mEdgeProperties.radiusVar = calculateVariance(mEdgeProperties.radii) / mDetectionVariables.predictedCircumference; // relative variance
            
            calculateCurvatureStats(mDetectionParameters, mEdgeProperties);
            
            mEdgeProperties.index     = iEdge;
            mEdgeProperties.tag       = 0;
            
            restoreEdgePoints(mEdgeProperties, cannyEdgesOriginal, cannyAOI); // Restore some points
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
    
    std::vector<edgeProperties> vEdgeCollectionProperties = edgeCollectionFilter(mDetectionVariables, mDetectionParameters, vEdgePropertiesNew, cannyAOI);

    std::vector<ellipseProperties> vEllipsePropertiesAll = ellipseFitting(mDetectionVariables, mDetectionParameters, vEdgeCollectionProperties, cannyAOI); // ellipse fitting
    
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
                
                if (std::isfinite(mEdgeProperties.curvature)) // ignore edges for which no curvature information is available
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
    
    mDetectionVariablesNew.thresholdChangeAspectRatio   = rangeChangeThresholdAspectRatio   * (1 - mDetectionVariables.certaintyFeatures) + mDetectionParameters.thresholdChangeAspectRatio;
    mDetectionVariablesNew.thresholdChangeCircumference = rangeChangeThresholdCircumference * (1 - mDetectionVariables.certaintyFeatures) + mDetectionParameters.thresholdChangeCircumference;
    mDetectionVariablesNew.thresholdChangePosition      = rangeChangeThresholdPosition      * (1 - mDetectionVariables.certaintyPosition) + mDetectionParameters.thresholdChangePosition;
    
    mDetectionVariablesNew.thresholdScore = mDetectionVariables.certaintyFeatures * mDetectionVariables.certaintyPosition * mDetectionParameters.thresholdScore;
    
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
    
    mDrawVariables.cannyEdgeIndices    = edgePointsOriginal;
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
        static const std::vector<double> arrayCurvatureMax =    {0.4125,0.3876,0.3650,0.3447,0.3266,0.3106,0.2965,0.2842,0.2734,0.2639,0.2556,0.2484,0.2419,0.2362,0.2310,0.2264,0.2222,0.2184,0.2148,0.2115,0.2084,0.2054,0.2026,0.1999,0.1973,0.1948,0.1923,0.1898,0.1874,0.1849,0.1824,0.1799,0.1774,0.1748,0.1722,0.1696,0.1669,0.1641,0.1614,0.1586,0.1559,0.1532,0.1506,0.1481,0.1456,0.1434,0.1413,0.1395,0.1378,0.1365,
                                                                 0.3984,0.3740,0.3518,0.3320,0.3144,0.2988,0.2852,0.2733,0.2628,0.2537,0.2458,0.2388,0.2326,0.2272,0.2223,0.2179,0.2140,0.2103,0.2070,0.2039,0.2010,0.1983,0.1957,0.1932,0.1907,0.1884,0.1860,0.1837,0.1814,0.1791,0.1768,0.1744,0.1720,0.1696,0.1671,0.1646,0.1620,0.1594,0.1568,0.1541,0.1515,0.1490,0.1465,0.1441,0.1418,0.1397,0.1378,0.1361,0.1346,0.1333,
                                                                 0.3849,0.3610,0.3393,0.3200,0.3029,0.2878,0.2746,0.2631,0.2530,0.2443,0.2367,0.2300,0.2242,0.2190,0.2144,0.2102,0.2065,0.2031,0.2000,0.1971,0.1944,0.1919,0.1895,0.1871,0.1849,0.1827,0.1805,0.1784,0.1762,0.1740,0.1718,0.1696,0.1673,0.1650,0.1626,0.1602,0.1577,0.1552,0.1527,0.1501,0.1476,0.1452,0.1428,0.1406,0.1384,0.1365,0.1347,0.1331,0.1317,0.1306,
                                                                 0.3721,0.3487,0.3276,0.3088,0.2921,0.2775,0.2648,0.2537,0.2440,0.2357,0.2284,0.2220,0.2164,0.2115,0.2072,0.2033,0.1998,0.1967,0.1938,0.1911,0.1886,0.1862,0.1840,0.1818,0.1797,0.1777,0.1757,0.1737,0.1716,0.1696,0.1675,0.1654,0.1632,0.1609,0.1586,0.1563,0.1539,0.1515,0.1490,0.1466,0.1442,0.1419,0.1396,0.1374,0.1354,0.1336,0.1319,0.1304,0.1292,0.1281,
                                                                 0.3601,0.3372,0.3166,0.2983,0.2822,0.2680,0.2557,0.2451,0.2358,0.2278,0.2208,0.2148,0.2095,0.2049,0.2008,0.1972,0.1939,0.1909,0.1882,0.1858,0.1834,0.1813,0.1792,0.1772,0.1752,0.1733,0.1714,0.1695,0.1676,0.1657,0.1637,0.1616,0.1595,0.1573,0.1551,0.1528,0.1505,0.1482,0.1458,0.1434,0.1411,0.1389,0.1367,0.1346,0.1327,0.1310,0.1294,0.1280,0.1269,0.1259,
                                                                 0.3489,0.3265,0.3065,0.2887,0.2730,0.2594,0.2475,0.2372,0.2284,0.2207,0.2141,0.2083,0.2033,0.1990,0.1951,0.1917,0.1887,0.1859,0.1834,0.1811,0.1790,0.1769,0.1750,0.1731,0.1713,0.1695,0.1677,0.1659,0.1641,0.1622,0.1603,0.1583,0.1563,0.1541,0.1520,0.1497,0.1475,0.1452,0.1429,0.1406,0.1383,0.1361,0.1341,0.1321,0.1302,0.1286,0.1271,0.1258,0.1247,0.1239,
                                                                 0.3385,0.3166,0.2971,0.2798,0.2647,0.2515,0.2401,0.2302,0.2217,0.2144,0.2081,0.2026,0.1979,0.1938,0.1901,0.1869,0.1841,0.1815,0.1792,0.1770,0.1750,0.1731,0.1713,0.1696,0.1679,0.1662,0.1644,0.1627,0.1610,0.1591,0.1573,0.1553,0.1533,0.1512,0.1491,0.1469,0.1447,0.1424,0.1402,0.1379,0.1357,0.1336,0.1316,0.1297,0.1280,0.1264,0.1250,0.1238,0.1228,0.1220,
                                                                 0.3289,0.3076,0.2886,0.2718,0.2571,0.2444,0.2334,0.2239,0.2158,0.2088,0.2028,0.1976,0.1931,0.1892,0.1858,0.1828,0.1801,0.1777,0.1755,0.1735,0.1716,0.1698,0.1681,0.1665,0.1648,0.1632,0.1615,0.1599,0.1582,0.1564,0.1545,0.1526,0.1506,0.1486,0.1465,0.1443,0.1421,0.1399,0.1377,0.1355,0.1333,0.1313,0.1293,0.1275,0.1258,0.1243,0.1230,0.1218,0.1209,0.1202,
                                                                 0.3201,0.2993,0.2808,0.2645,0.2503,0.2380,0.2274,0.2183,0.2105,0.2038,0.1980,0.1931,0.1889,0.1852,0.1819,0.1791,0.1766,0.1743,0.1723,0.1704,0.1686,0.1669,0.1653,0.1637,0.1621,0.1605,0.1589,0.1573,0.1556,0.1539,0.1520,0.1501,0.1482,0.1461,0.1440,0.1419,0.1397,0.1375,0.1353,0.1331,0.1311,0.1291,0.1272,0.1254,0.1238,0.1223,0.1211,0.1200,0.1191,0.1184,
                                                                 0.3120,0.2917,0.2737,0.2579,0.2442,0.2323,0.2221,0.2133,0.2058,0.1994,0.1939,0.1892,0.1851,0.1816,0.1786,0.1759,0.1735,0.1713,0.1694,0.1676,0.1659,0.1643,0.1627,0.1612,0.1597,0.1581,0.1565,0.1549,0.1532,0.1515,0.1497,0.1478,0.1458,0.1438,0.1417,0.1395,0.1374,0.1352,0.1330,0.1309,0.1289,0.1269,0.1251,0.1234,0.1218,0.1205,0.1193,0.1182,0.1174,0.1168,
                                                                 0.3046,0.2848,0.2673,0.2520,0.2386,0.2271,0.2173,0.2088,0.2016,0.1954,0.1902,0.1857,0.1818,0.1785,0.1756,0.1730,0.1707,0.1687,0.1668,0.1651,0.1635,0.1619,0.1604,0.1589,0.1574,0.1559,0.1543,0.1527,0.1510,0.1493,0.1475,0.1456,0.1436,0.1415,0.1394,0.1373,0.1351,0.1330,0.1309,0.1288,0.1268,0.1249,0.1231,0.1215,0.1200,0.1187,0.1175,0.1166,0.1158,0.1152,
                                                                 0.2978,0.2785,0.2614,0.2465,0.2336,0.2225,0.2129,0.2048,0.1978,0.1919,0.1869,0.1826,0.1789,0.1757,0.1729,0.1704,0.1682,0.1663,0.1645,0.1628,0.1612,0.1597,0.1582,0.1568,0.1553,0.1538,0.1522,0.1506,0.1489,0.1472,0.1453,0.1434,0.1414,0.1394,0.1373,0.1351,0.1330,0.1308,0.1287,0.1267,0.1248,0.1229,0.1212,0.1196,0.1182,0.1169,0.1158,0.1149,0.1142,0.1137,
                                                                 0.2915,0.2727,0.2561,0.2416,0.2290,0.2182,0.2090,0.2011,0.1944,0.1887,0.1839,0.1797,0.1762,0.1731,0.1704,0.1681,0.1660,0.1641,0.1623,0.1607,0.1592,0.1577,0.1562,0.1548,0.1533,0.1518,0.1502,0.1486,0.1469,0.1451,0.1432,0.1413,0.1393,0.1372,0.1351,0.1330,0.1309,0.1288,0.1267,0.1247,0.1228,0.1210,0.1193,0.1178,0.1165,0.1153,0.1142,0.1134,0.1127,0.1123,
                                                                 0.2857,0.2673,0.2511,0.2370,0.2248,0.2143,0.2054,0.1978,0.1913,0.1858,0.1811,0.1771,0.1737,0.1707,0.1681,0.1659,0.1638,0.1620,0.1603,0.1587,0.1572,0.1557,0.1543,0.1528,0.1514,0.1498,0.1483,0.1466,0.1449,0.1431,0.1412,0.1393,0.1372,0.1352,0.1331,0.1309,0.1288,0.1267,0.1247,0.1228,0.1209,0.1192,0.1176,0.1161,0.1148,0.1137,0.1127,0.1119,0.1113,0.1109,
                                                                 0.2803,0.2623,0.2465,0.2327,0.2209,0.2107,0.2021,0.1947,0.1884,0.1831,0.1786,0.1747,0.1714,0.1685,0.1660,0.1638,0.1618,0.1600,0.1584,0.1568,0.1553,0.1539,0.1524,0.1510,0.1495,0.1480,0.1464,0.1447,0.1430,0.1411,0.1392,0.1373,0.1352,0.1331,0.1310,0.1289,0.1268,0.1248,0.1228,0.1209,0.1191,0.1174,0.1159,0.1145,0.1132,0.1121,0.1112,0.1105,0.1099,0.1095,
                                                                 0.2752,0.2576,0.2421,0.2288,0.2172,0.2074,0.1989,0.1918,0.1857,0.1806,0.1762,0.1724,0.1692,0.1664,0.1640,0.1619,0.1599,0.1582,0.1566,0.1550,0.1535,0.1521,0.1506,0.1492,0.1477,0.1461,0.1445,0.1428,0.1410,0.1392,0.1373,0.1353,0.1332,0.1312,0.1291,0.1270,0.1249,0.1229,0.1209,0.1191,0.1173,0.1157,0.1142,0.1129,0.1117,0.1107,0.1098,0.1091,0.1086,0.1083,
                                                                 0.2703,0.2531,0.2380,0.2250,0.2138,0.2042,0.1960,0.1891,0.1832,0.1782,0.1739,0.1703,0.1672,0.1645,0.1621,0.1600,0.1581,0.1564,0.1548,0.1533,0.1518,0.1503,0.1489,0.1474,0.1459,0.1443,0.1427,0.1410,0.1392,0.1373,0.1354,0.1333,0.1313,0.1292,0.1271,0.1250,0.1230,0.1210,0.1191,0.1173,0.1156,0.1141,0.1126,0.1114,0.1103,0.1093,0.1085,0.1079,0.1074,0.1071,
                                                                 0.2657,0.2488,0.2341,0.2214,0.2105,0.2011,0.1932,0.1864,0.1807,0.1759,0.1717,0.1682,0.1652,0.1626,0.1603,0.1582,0.1564,0.1547,0.1531,0.1516,0.1501,0.1487,0.1472,0.1457,0.1442,0.1426,0.1409,0.1391,0.1373,0.1354,0.1335,0.1315,0.1294,0.1273,0.1252,0.1232,0.1212,0.1192,0.1174,0.1156,0.1140,0.1125,0.1112,0.1099,0.1089,0.1080,0.1072,0.1066,0.1062,0.1059,
                                                                 0.2612,0.2447,0.2304,0.2180,0.2073,0.1982,0.1905,0.1840,0.1784,0.1737,0.1697,0.1662,0.1633,0.1607,0.1585,0.1565,0.1547,0.1530,0.1514,0.1499,0.1485,0.1470,0.1455,0.1440,0.1425,0.1408,0.1391,0.1374,0.1355,0.1336,0.1316,0.1296,0.1275,0.1255,0.1234,0.1214,0.1194,0.1175,0.1157,0.1140,0.1125,0.1110,0.1097,0.1086,0.1076,0.1067,0.1060,0.1055,0.1051,0.1049,
                                                                 0.2569,0.2408,0.2268,0.2147,0.2043,0.1954,0.1879,0.1816,0.1762,0.1716,0.1677,0.1643,0.1615,0.1590,0.1568,0.1548,0.1530,0.1514,0.1498,0.1483,0.1469,0.1454,0.1439,0.1424,0.1408,0.1391,0.1374,0.1356,0.1337,0.1318,0.1298,0.1278,0.1257,0.1237,0.1216,0.1196,0.1177,0.1159,0.1141,0.1125,0.1110,0.1096,0.1084,0.1073,0.1064,0.1056,0.1049,0.1044,0.1041,0.1039,
                                                                 0.2527,0.2370,0.2233,0.2115,0.2014,0.1928,0.1855,0.1793,0.1740,0.1696,0.1658,0.1625,0.1597,0.1573,0.1551,0.1532,0.1514,0.1498,0.1483,0.1468,0.1453,0.1438,0.1423,0.1408,0.1391,0.1375,0.1357,0.1339,0.1320,0.1301,0.1281,0.1260,0.1240,0.1219,0.1199,0.1180,0.1161,0.1143,0.1126,0.1110,0.1096,0.1083,0.1071,0.1061,0.1052,0.1045,0.1039,0.1034,0.1031,0.1029,
                                                                 0.2487,0.2333,0.2199,0.2084,0.1986,0.1902,0.1831,0.1770,0.1719,0.1676,0.1639,0.1608,0.1580,0.1557,0.1535,0.1517,0.1499,0.1483,0.1468,0.1453,0.1438,0.1423,0.1408,0.1392,0.1376,0.1359,0.1341,0.1322,0.1303,0.1284,0.1264,0.1243,0.1223,0.1203,0.1183,0.1164,0.1145,0.1128,0.1112,0.1097,0.1083,0.1070,0.1059,0.1050,0.1041,0.1035,0.1029,0.1025,0.1022,0.1021,
                                                                 0.2448,0.2297,0.2167,0.2054,0.1958,0.1877,0.1808,0.1749,0.1700,0.1657,0.1622,0.1591,0.1564,0.1541,0.1520,0.1502,0.1484,0.1468,0.1453,0.1438,0.1423,0.1408,0.1392,0.1377,0.1360,0.1343,0.1325,0.1306,0.1287,0.1267,0.1247,0.1227,0.1206,0.1186,0.1167,0.1148,0.1130,0.1114,0.1098,0.1084,0.1070,0.1059,0.1048,0.1039,0.1031,0.1025,0.1020,0.1016,0.1014,0.1013,
                                                                 0.2410,0.2263,0.2135,0.2026,0.1932,0.1853,0.1786,0.1729,0.1680,0.1640,0.1605,0.1575,0.1549,0.1526,0.1506,0.1487,0.1470,0.1454,0.1439,0.1424,0.1409,0.1393,0.1378,0.1362,0.1345,0.1327,0.1309,0.1290,0.1271,0.1251,0.1231,0.1211,0.1191,0.1171,0.1152,0.1134,0.1116,0.1100,0.1085,0.1071,0.1059,0.1048,0.1038,0.1029,0.1022,0.1016,0.1011,0.1008,0.1006,0.1005,
                                                                 0.2374,0.2230,0.2105,0.1998,0.1907,0.1830,0.1765,0.1709,0.1662,0.1622,0.1588,0.1559,0.1534,0.1511,0.1491,0.1473,0.1457,0.1441,0.1425,0.1410,0.1395,0.1379,0.1364,0.1347,0.1330,0.1312,0.1294,0.1275,0.1255,0.1235,0.1215,0.1195,0.1175,0.1156,0.1137,0.1120,0.1103,0.1087,0.1073,0.1060,0.1048,0.1037,0.1028,0.1020,0.1013,0.1008,0.1003,0.1000,0.0999,0.0998,
                                                                 0.2338,0.2197,0.2076,0.1972,0.1883,0.1808,0.1744,0.1690,0.1645,0.1606,0.1573,0.1544,0.1519,0.1498,0.1478,0.1460,0.1443,0.1427,0.1412,0.1397,0.1381,0.1366,0.1350,0.1333,0.1315,0.1297,0.1279,0.1260,0.1240,0.1220,0.1200,0.1180,0.1161,0.1142,0.1124,0.1106,0.1090,0.1075,0.1061,0.1049,0.1038,0.1028,0.1019,0.1011,0.1005,0.1000,0.0996,0.0993,0.0992,0.0992,
                                                                 0.2304,0.2166,0.2048,0.1946,0.1860,0.1787,0.1725,0.1673,0.1628,0.1590,0.1558,0.1530,0.1506,0.1484,0.1465,0.1447,0.1431,0.1415,0.1399,0.1384,0.1368,0.1352,0.1336,0.1319,0.1301,0.1283,0.1264,0.1245,0.1225,0.1205,0.1185,0.1166,0.1147,0.1128,0.1111,0.1094,0.1078,0.1064,0.1051,0.1039,0.1028,0.1019,0.1010,0.1003,0.0998,0.0993,0.0989,0.0987,0.0986,0.0986,
                                                                 0.2271,0.2137,0.2021,0.1922,0.1838,0.1767,0.1706,0.1655,0.1612,0.1575,0.1544,0.1517,0.1493,0.1472,0.1453,0.1435,0.1419,0.1403,0.1387,0.1371,0.1356,0.1340,0.1323,0.1306,0.1288,0.1269,0.1250,0.1231,0.1211,0.1191,0.1171,0.1152,0.1133,0.1115,0.1098,0.1082,0.1067,0.1053,0.1040,0.1029,0.1019,0.1010,0.1003,0.0996,0.0991,0.0986,0.0983,0.0981,0.0980,0.0980,
                                                                 0.2239,0.2108,0.1995,0.1899,0.1817,0.1747,0.1689,0.1639,0.1597,0.1561,0.1530,0.1504,0.1480,0.1460,0.1441,0.1423,0.1407,0.1391,0.1375,0.1359,0.1343,0.1327,0.1310,0.1292,0.1274,0.1256,0.1236,0.1217,0.1197,0.1177,0.1158,0.1139,0.1120,0.1103,0.1086,0.1071,0.1056,0.1043,0.1031,0.1020,0.1011,0.1002,0.0995,0.0989,0.0984,0.0980,0.0977,0.0975,0.0975,0.0975,
                                                                 0.2209,0.2080,0.1970,0.1876,0.1797,0.1729,0.1672,0.1624,0.1583,0.1548,0.1518,0.1492,0.1469,0.1448,0.1429,0.1412,0.1395,0.1379,0.1364,0.1348,0.1331,0.1315,0.1297,0.1280,0.1261,0.1242,0.1223,0.1203,0.1184,0.1164,0.1145,0.1126,0.1108,0.1091,0.1075,0.1060,0.1046,0.1034,0.1022,0.1012,0.1003,0.0995,0.0988,0.0983,0.0978,0.0975,0.0972,0.0970,0.0970,0.0970,
                                                                 0.2179,0.2054,0.1946,0.1855,0.1777,0.1712,0.1656,0.1609,0.1569,0.1535,0.1505,0.1480,0.1457,0.1437,0.1418,0.1401,0.1385,0.1368,0.1352,0.1336,0.1320,0.1303,0.1285,0.1267,0.1248,0.1229,0.1210,0.1190,0.1171,0.1151,0.1132,0.1114,0.1096,0.1080,0.1064,0.1050,0.1037,0.1025,0.1014,0.1004,0.0996,0.0988,0.0982,0.0977,0.0973,0.0969,0.0967,0.0966,0.0965,0.0966,
                                                                 0.2151,0.2028,0.1924,0.1834,0.1759,0.1695,0.1641,0.1595,0.1556,0.1522,0.1494,0.1469,0.1446,0.1426,0.1408,0.1391,0.1374,0.1358,0.1341,0.1325,0.1308,0.1291,0.1273,0.1255,0.1236,0.1217,0.1197,0.1178,0.1158,0.1139,0.1120,0.1102,0.1085,0.1069,0.1054,0.1040,0.1028,0.1016,0.1006,0.0997,0.0989,0.0982,0.0976,0.0971,0.0968,0.0965,0.0962,0.0961,0.0961,0.0962,
                                                                 0.2123,0.2004,0.1902,0.1815,0.1741,0.1679,0.1626,0.1581,0.1543,0.1511,0.1483,0.1458,0.1436,0.1416,0.1398,0.1380,0.1364,0.1347,0.1331,0.1314,0.1297,0.1279,0.1261,0.1243,0.1224,0.1204,0.1185,0.1165,0.1146,0.1127,0.1109,0.1091,0.1075,0.1059,0.1045,0.1031,0.1019,0.1009,0.0999,0.0990,0.0983,0.0976,0.0971,0.0966,0.0963,0.0960,0.0958,0.0957,0.0957,0.0958,
                                                                 0.2097,0.1980,0.1881,0.1796,0.1724,0.1664,0.1612,0.1569,0.1531,0.1500,0.1472,0.1448,0.1426,0.1406,0.1388,0.1370,0.1354,0.1337,0.1320,0.1303,0.1286,0.1268,0.1250,0.1231,0.1212,0.1192,0.1173,0.1153,0.1134,0.1116,0.1098,0.1081,0.1065,0.1050,0.1036,0.1023,0.1012,0.1001,0.0992,0.0984,0.0977,0.0971,0.0966,0.0962,0.0958,0.0956,0.0954,0.0953,0.0953,0.0954,
                                                                 0.2072,0.1958,0.1861,0.1778,0.1708,0.1649,0.1599,0.1556,0.1520,0.1489,0.1462,0.1438,0.1416,0.1397,0.1378,0.1361,0.1344,0.1327,0.1310,0.1293,0.1275,0.1257,0.1238,0.1219,0.1200,0.1180,0.1161,0.1142,0.1123,0.1105,0.1087,0.1070,0.1055,0.1041,0.1027,0.1015,0.1004,0.0995,0.0986,0.0978,0.0972,0.0966,0.0961,0.0958,0.0954,0.0952,0.0951,0.0950,0.0950,0.0951,
                                                                 0.2047,0.1936,0.1841,0.1761,0.1693,0.1635,0.1586,0.1545,0.1509,0.1479,0.1452,0.1428,0.1407,0.1387,0.1369,0.1351,0.1334,0.1317,0.1300,0.1282,0.1264,0.1246,0.1227,0.1208,0.1188,0.1169,0.1150,0.1130,0.1112,0.1094,0.1077,0.1061,0.1046,0.1032,0.1019,0.1008,0.0997,0.0988,0.0980,0.0973,0.0967,0.0962,0.0957,0.0954,0.0951,0.0949,0.0947,0.0947,0.0947,0.0948,
                                                                 0.2024,0.1915,0.1823,0.1744,0.1678,0.1622,0.1574,0.1533,0.1499,0.1469,0.1442,0.1419,0.1398,0.1378,0.1360,0.1342,0.1325,0.1307,0.1290,0.1272,0.1254,0.1235,0.1216,0.1197,0.1177,0.1158,0.1138,0.1120,0.1101,0.1084,0.1067,0.1052,0.1037,0.1024,0.1012,0.1001,0.0991,0.0982,0.0975,0.0968,0.0962,0.0957,0.0953,0.0950,0.0947,0.0946,0.0944,0.0944,0.0944,0.0946,
                                                                 0.2001,0.1895,0.1805,0.1728,0.1664,0.1609,0.1562,0.1523,0.1489,0.1459,0.1433,0.1410,0.1389,0.1370,0.1351,0.1333,0.1315,0.1298,0.1280,0.1262,0.1243,0.1224,0.1205,0.1186,0.1166,0.1147,0.1128,0.1109,0.1091,0.1074,0.1058,0.1043,0.1029,0.1016,0.1005,0.0994,0.0985,0.0977,0.0970,0.0963,0.0958,0.0953,0.0950,0.0947,0.0944,0.0943,0.0942,0.0941,0.0942,0.0943,
                                                                 0.1979,0.1876,0.1788,0.1713,0.1650,0.1597,0.1551,0.1512,0.1479,0.1450,0.1424,0.1401,0.1380,0.1361,0.1342,0.1324,0.1306,0.1288,0.1270,0.1252,0.1233,0.1214,0.1194,0.1175,0.1155,0.1136,0.1117,0.1099,0.1081,0.1065,0.1049,0.1035,0.1021,0.1009,0.0998,0.0988,0.0979,0.0972,0.0965,0.0959,0.0954,0.0950,0.0946,0.0944,0.0941,0.0940,0.0939,0.0939,0.0939,0.0941,
                                                                 0.1958,0.1857,0.1771,0.1699,0.1637,0.1585,0.1540,0.1502,0.1470,0.1441,0.1416,0.1393,0.1372,0.1352,0.1334,0.1315,0.1297,0.1279,0.1260,0.1242,0.1223,0.1203,0.1184,0.1164,0.1145,0.1126,0.1107,0.1089,0.1072,0.1056,0.1041,0.1027,0.1014,0.1002,0.0992,0.0983,0.0974,0.0967,0.0961,0.0955,0.0950,0.0947,0.0943,0.0941,0.0939,0.0937,0.0937,0.0937,0.0937,0.0939,
                                                                 0.1938,0.1839,0.1756,0.1685,0.1624,0.1573,0.1530,0.1493,0.1461,0.1432,0.1407,0.1385,0.1364,0.1344,0.1325,0.1306,0.1288,0.1270,0.1251,0.1232,0.1213,0.1193,0.1173,0.1154,0.1135,0.1116,0.1097,0.1080,0.1063,0.1047,0.1033,0.1019,0.1007,0.0996,0.0986,0.0977,0.0969,0.0963,0.0957,0.0951,0.0947,0.0943,0.0940,0.0938,0.0936,0.0935,0.0935,0.0935,0.0935,0.0937,
                                                                 0.1918,0.1822,0.1740,0.1671,0.1612,0.1562,0.1520,0.1483,0.1452,0.1424,0.1399,0.1376,0.1356,0.1336,0.1317,0.1298,0.1279,0.1260,0.1241,0.1222,0.1202,0.1183,0.1163,0.1144,0.1124,0.1106,0.1088,0.1071,0.1054,0.1039,0.1025,0.1012,0.1001,0.0990,0.0981,0.0972,0.0965,0.0958,0.0953,0.0948,0.0944,0.0941,0.0938,0.0936,0.0934,0.0933,0.0933,0.0933,0.0933,0.0935,
                                                                 0.1899,0.1805,0.1725,0.1658,0.1600,0.1552,0.1510,0.1474,0.1443,0.1416,0.1391,0.1368,0.1347,0.1327,0.1308,0.1289,0.1270,0.1251,0.1232,0.1212,0.1192,0.1173,0.1153,0.1134,0.1115,0.1096,0.1079,0.1062,0.1046,0.1031,0.1018,0.1006,0.0994,0.0984,0.0976,0.0968,0.0961,0.0955,0.0949,0.0945,0.0941,0.0938,0.0935,0.0933,0.0932,0.0931,0.0931,0.0931,0.0932,0.0933,
                                                                 0.1881,0.1789,0.1711,0.1645,0.1589,0.1541,0.1501,0.1465,0.1435,0.1407,0.1383,0.1360,0.1339,0.1319,0.1300,0.1280,0.1261,0.1242,0.1222,0.1202,0.1183,0.1163,0.1143,0.1124,0.1105,0.1087,0.1070,0.1053,0.1038,0.1024,0.1011,0.0999,0.0989,0.0979,0.0971,0.0963,0.0957,0.0951,0.0946,0.0942,0.0938,0.0936,0.0933,0.0931,0.0930,0.0929,0.0929,0.0929,0.0930,0.0932,
                                                                 0.1863,0.1773,0.1697,0.1633,0.1578,0.1531,0.1491,0.1457,0.1426,0.1399,0.1375,0.1352,0.1331,0.1311,0.1291,0.1271,0.1252,0.1232,0.1212,0.1193,0.1173,0.1153,0.1133,0.1114,0.1096,0.1078,0.1061,0.1045,0.1031,0.1017,0.1005,0.0993,0.0983,0.0974,0.0966,0.0959,0.0953,0.0948,0.0943,0.0939,0.0936,0.0933,0.0931,0.0930,0.0928,0.0928,0.0927,0.0928,0.0929,0.0930,
                                                                 0.1845,0.1758,0.1684,0.1621,0.1567,0.1522,0.1482,0.1448,0.1418,0.1391,0.1367,0.1344,0.1323,0.1303,0.1283,0.1263,0.1243,0.1223,0.1203,0.1183,0.1163,0.1143,0.1124,0.1105,0.1087,0.1069,0.1053,0.1038,0.1023,0.1010,0.0998,0.0988,0.0978,0.0970,0.0962,0.0955,0.0950,0.0945,0.0940,0.0937,0.0934,0.0931,0.0929,0.0928,0.0927,0.0926,0.0926,0.0926,0.0927,0.0929,
                                                                 0.1829,0.1743,0.1671,0.1609,0.1557,0.1512,0.1473,0.1440,0.1410,0.1383,0.1359,0.1336,0.1315,0.1294,0.1274,0.1254,0.1234,0.1214,0.1193,0.1173,0.1153,0.1134,0.1114,0.1096,0.1078,0.1061,0.1045,0.1030,0.1016,0.1004,0.0993,0.0982,0.0973,0.0965,0.0958,0.0952,0.0947,0.0942,0.0938,0.0935,0.0932,0.0929,0.0928,0.0926,0.0925,0.0925,0.0925,0.0925,0.0926,0.0927,
                                                                 0.1812,0.1729,0.1658,0.1598,0.1547,0.1503,0.1465,0.1431,0.1402,0.1376,0.1351,0.1329,0.1307,0.1286,0.1265,0.1245,0.1225,0.1204,0.1184,0.1164,0.1144,0.1124,0.1105,0.1087,0.1069,0.1053,0.1037,0.1023,0.1010,0.0998,0.0987,0.0977,0.0969,0.0961,0.0955,0.0949,0.0944,0.0939,0.0936,0.0932,0.0930,0.0928,0.0926,0.0925,0.0924,0.0923,0.0923,0.0924,0.0925,0.0926,
                                                                 0.1796,0.1715,0.1646,0.1587,0.1537,0.1493,0.1456,0.1423,0.1394,0.1368,0.1343,0.1320,0.1299,0.1277,0.1257,0.1236,0.1215,0.1195,0.1175,0.1154,0.1134,0.1115,0.1096,0.1078,0.1061,0.1045,0.1030,0.1016,0.1004,0.0992,0.0982,0.0973,0.0965,0.0957,0.0951,0.0946,0.0941,0.0937,0.0933,0.0930,0.0928,0.0926,0.0924,0.0923,0.0923,0.0922,0.0922,0.0923,0.0924,0.0925,
                                                                 0.1781,0.1701,0.1634,0.1576,0.1527,0.1484,0.1447,0.1415,0.1386,0.1360,0.1335,0.1312,0.1290,0.1269,0.1248,0.1227,0.1206,0.1186,0.1165,0.1145,0.1125,0.1106,0.1087,0.1070,0.1053,0.1037,0.1023,0.1010,0.0998,0.0987,0.0977,0.0968,0.0961,0.0954,0.0948,0.0943,0.0938,0.0935,0.0931,0.0929,0.0926,0.0925,0.0923,0.0922,0.0921,0.0921,0.0921,0.0922,0.0923,0.0924};
        
        return arrayCurvatureMax[y * arrayWidth + x];
    }
    case 6:
    {
        static const std::vector<double> arrayCurvatureMax =    {0.3533,0.3385,0.3240,0.3100,0.2965,0.2838,0.2718,0.2606,0.2502,0.2406,0.2317,0.2235,0.2160,0.2091,0.2028,0.1970,0.1917,0.1869,0.1825,0.1785,0.1748,0.1715,0.1685,0.1657,0.1632,0.1609,0.1587,0.1567,0.1548,0.1529,0.1512,0.1495,0.1478,0.1463,0.1447,0.1433,0.1420,0.1408,0.1398,0.1389,0.1381,0.1374,0.1368,0.1364,0.1360,0.1356,0.1353,0.1350,0.1348,0.1345,
                                                                 0.3470,0.3321,0.3174,0.3032,0.2896,0.2767,0.2647,0.2534,0.2430,0.2334,0.2246,0.2165,0.2091,0.2023,0.1961,0.1905,0.1853,0.1807,0.1765,0.1726,0.1692,0.1661,0.1633,0.1607,0.1584,0.1562,0.1543,0.1524,0.1507,0.1490,0.1474,0.1459,0.1444,0.1430,0.1417,0.1404,0.1393,0.1383,0.1374,0.1366,0.1359,0.1354,0.1349,0.1346,0.1343,0.1340,0.1338,0.1336,0.1334,0.1332,
                                                                 0.3406,0.3255,0.3106,0.2963,0.2826,0.2697,0.2576,0.2463,0.2360,0.2264,0.2177,0.2097,0.2024,0.1958,0.1898,0.1843,0.1794,0.1749,0.1709,0.1673,0.1640,0.1611,0.1584,0.1561,0.1539,0.1519,0.1501,0.1484,0.1468,0.1453,0.1439,0.1425,0.1411,0.1398,0.1386,0.1375,0.1365,0.1356,0.1348,0.1342,0.1337,0.1332,0.1329,0.1326,0.1324,0.1323,0.1321,0.1320,0.1318,0.1317,
                                                                 0.3340,0.3187,0.3037,0.2893,0.2756,0.2626,0.2506,0.2394,0.2291,0.2197,0.2111,0.2032,0.1961,0.1897,0.1839,0.1786,0.1739,0.1696,0.1658,0.1624,0.1593,0.1565,0.1541,0.1518,0.1498,0.1480,0.1463,0.1447,0.1432,0.1418,0.1405,0.1392,0.1379,0.1368,0.1356,0.1346,0.1337,0.1329,0.1323,0.1317,0.1313,0.1310,0.1307,0.1306,0.1304,0.1303,0.1303,0.1302,0.1301,0.1300,
                                                                 0.3273,0.3118,0.2967,0.2823,0.2685,0.2557,0.2437,0.2326,0.2225,0.2132,0.2048,0.1971,0.1902,0.1840,0.1784,0.1733,0.1688,0.1648,0.1611,0.1579,0.1550,0.1524,0.1501,0.1480,0.1461,0.1444,0.1428,0.1413,0.1399,0.1386,0.1373,0.1361,0.1349,0.1338,0.1327,0.1318,0.1310,0.1303,0.1297,0.1292,0.1289,0.1286,0.1284,0.1283,0.1283,0.1282,0.1282,0.1282,0.1281,0.1280,
                                                                 0.3204,0.3048,0.2897,0.2752,0.2616,0.2488,0.2369,0.2260,0.2160,0.2070,0.1988,0.1913,0.1847,0.1787,0.1733,0.1685,0.1642,0.1604,0.1569,0.1539,0.1511,0.1487,0.1465,0.1446,0.1428,0.1412,0.1397,0.1383,0.1369,0.1356,0.1344,0.1332,0.1320,0.1310,0.1300,0.1291,0.1283,0.1276,0.1271,0.1266,0.1263,0.1261,0.1260,0.1260,0.1260,0.1260,0.1260,0.1260,0.1260,0.1259,
                                                                 0.3134,0.2978,0.2827,0.2683,0.2547,0.2420,0.2303,0.2196,0.2099,0.2011,0.1931,0.1859,0.1795,0.1738,0.1686,0.1641,0.1600,0.1564,0.1531,0.1503,0.1477,0.1454,0.1434,0.1415,0.1398,0.1383,0.1368,0.1355,0.1342,0.1329,0.1317,0.1305,0.1294,0.1283,0.1273,0.1264,0.1257,0.1250,0.1245,0.1241,0.1238,0.1236,0.1235,0.1235,0.1235,0.1235,0.1236,0.1236,0.1236,0.1235,
                                                                 0.3064,0.2907,0.2757,0.2614,0.2479,0.2355,0.2240,0.2135,0.2040,0.1955,0.1878,0.1809,0.1747,0.1692,0.1644,0.1600,0.1562,0.1528,0.1497,0.1470,0.1446,0.1425,0.1406,0.1388,0.1372,0.1357,0.1343,0.1330,0.1317,0.1305,0.1293,0.1281,0.1270,0.1259,0.1249,0.1240,0.1232,0.1225,0.1220,0.1216,0.1213,0.1211,0.1210,0.1210,0.1210,0.1210,0.1210,0.1210,0.1210,0.1210,
                                                                 0.2993,0.2837,0.2688,0.2546,0.2414,0.2291,0.2179,0.2077,0.1985,0.1902,0.1828,0.1762,0.1703,0.1651,0.1605,0.1564,0.1528,0.1495,0.1467,0.1442,0.1419,0.1399,0.1381,0.1364,0.1349,0.1335,0.1321,0.1308,0.1296,0.1283,0.1271,0.1259,0.1247,0.1236,0.1226,0.1217,0.1208,0.1201,0.1196,0.1191,0.1188,0.1186,0.1184,0.1184,0.1183,0.1183,0.1183,0.1183,0.1183,0.1182,
                                                                 0.2923,0.2768,0.2620,0.2480,0.2350,0.2230,0.2121,0.2022,0.1932,0.1853,0.1781,0.1718,0.1662,0.1613,0.1569,0.1531,0.1497,0.1467,0.1440,0.1416,0.1395,0.1376,0.1359,0.1343,0.1329,0.1315,0.1302,0.1289,0.1277,0.1264,0.1252,0.1239,0.1227,0.1216,0.1205,0.1195,0.1187,0.1179,0.1173,0.1168,0.1164,0.1161,0.1159,0.1158,0.1157,0.1156,0.1155,0.1155,0.1154,0.1153,
                                                                 0.2853,0.2700,0.2553,0.2416,0.2289,0.2172,0.2065,0.1969,0.1883,0.1807,0.1739,0.1678,0.1625,0.1578,0.1537,0.1501,0.1469,0.1441,0.1416,0.1394,0.1374,0.1356,0.1340,0.1325,0.1311,0.1298,0.1285,0.1272,0.1260,0.1247,0.1234,0.1222,0.1209,0.1198,0.1186,0.1176,0.1167,0.1158,0.1151,0.1146,0.1141,0.1137,0.1134,0.1132,0.1130,0.1128,0.1127,0.1126,0.1124,0.1123,
                                                                 0.2785,0.2633,0.2489,0.2355,0.2230,0.2116,0.2013,0.1920,0.1838,0.1764,0.1699,0.1642,0.1591,0.1547,0.1508,0.1474,0.1444,0.1418,0.1394,0.1374,0.1355,0.1338,0.1323,0.1309,0.1295,0.1282,0.1270,0.1257,0.1245,0.1232,0.1219,0.1206,0.1193,0.1181,0.1169,0.1158,0.1148,0.1139,0.1131,0.1124,0.1119,0.1114,0.1110,0.1107,0.1104,0.1101,0.1099,0.1096,0.1094,0.1092,
                                                                 0.2718,0.2569,0.2427,0.2296,0.2174,0.2064,0.1964,0.1875,0.1795,0.1725,0.1663,0.1608,0.1561,0.1519,0.1482,0.1450,0.1422,0.1397,0.1375,0.1356,0.1338,0.1322,0.1308,0.1294,0.1281,0.1269,0.1256,0.1244,0.1231,0.1218,0.1205,0.1192,0.1179,0.1166,0.1153,0.1142,0.1131,0.1121,0.1112,0.1105,0.1098,0.1092,0.1087,0.1082,0.1078,0.1074,0.1071,0.1067,0.1064,0.1061,
                                                                 0.2653,0.2507,0.2368,0.2240,0.2122,0.2015,0.1918,0.1832,0.1756,0.1689,0.1629,0.1578,0.1533,0.1493,0.1459,0.1428,0.1402,0.1379,0.1358,0.1340,0.1324,0.1309,0.1295,0.1282,0.1269,0.1257,0.1245,0.1232,0.1219,0.1206,0.1193,0.1179,0.1165,0.1152,0.1139,0.1126,0.1115,0.1104,0.1094,0.1086,0.1078,0.1071,0.1064,0.1059,0.1053,0.1048,0.1044,0.1039,0.1035,0.1031,
                                                                 0.2591,0.2447,0.2312,0.2187,0.2072,0.1968,0.1875,0.1793,0.1720,0.1655,0.1599,0.1550,0.1507,0.1470,0.1437,0.1409,0.1384,0.1362,0.1343,0.1326,0.1310,0.1296,0.1283,0.1270,0.1258,0.1246,0.1234,0.1221,0.1209,0.1195,0.1181,0.1167,0.1153,0.1139,0.1125,0.1112,0.1100,0.1088,0.1078,0.1068,0.1059,0.1051,0.1043,0.1036,0.1029,0.1023,0.1017,0.1012,0.1007,0.1003,
                                                                 0.2532,0.2391,0.2259,0.2137,0.2026,0.1925,0.1836,0.1756,0.1686,0.1625,0.1571,0.1525,0.1484,0.1449,0.1418,0.1392,0.1368,0.1348,0.1330,0.1313,0.1299,0.1285,0.1272,0.1260,0.1248,0.1236,0.1224,0.1212,0.1199,0.1185,0.1171,0.1156,0.1142,0.1127,0.1113,0.1099,0.1086,0.1073,0.1061,0.1051,0.1041,0.1031,0.1022,0.1014,0.1006,0.0999,0.0992,0.0986,0.0980,0.0975,
                                                                 0.2475,0.2337,0.2208,0.2090,0.1982,0.1886,0.1799,0.1723,0.1656,0.1597,0.1546,0.1502,0.1463,0.1430,0.1401,0.1376,0.1354,0.1335,0.1318,0.1302,0.1288,0.1275,0.1263,0.1251,0.1239,0.1228,0.1215,0.1203,0.1189,0.1175,0.1161,0.1146,0.1131,0.1115,0.1100,0.1086,0.1072,0.1058,0.1046,0.1034,0.1023,0.1013,0.1003,0.0993,0.0985,0.0976,0.0969,0.0962,0.0956,0.0950,
                                                                 0.2421,0.2287,0.2161,0.2046,0.1942,0.1849,0.1765,0.1692,0.1628,0.1572,0.1523,0.1481,0.1445,0.1413,0.1386,0.1362,0.1341,0.1323,0.1307,0.1292,0.1279,0.1266,0.1254,0.1243,0.1231,0.1219,0.1207,0.1194,0.1181,0.1166,0.1151,0.1136,0.1120,0.1104,0.1088,0.1073,0.1058,0.1044,0.1031,0.1018,0.1006,0.0994,0.0984,0.0973,0.0964,0.0955,0.0947,0.0939,0.0933,0.0927,
                                                                 0.2371,0.2239,0.2117,0.2006,0.1905,0.1815,0.1734,0.1664,0.1602,0.1549,0.1502,0.1462,0.1428,0.1398,0.1372,0.1349,0.1330,0.1312,0.1297,0.1283,0.1270,0.1258,0.1246,0.1235,0.1224,0.1212,0.1199,0.1186,0.1172,0.1157,0.1142,0.1126,0.1109,0.1093,0.1076,0.1060,0.1045,0.1030,0.1015,0.1002,0.0989,0.0977,0.0965,0.0955,0.0944,0.0935,0.0926,0.0919,0.0912,0.0905,
                                                                 0.2323,0.2195,0.2076,0.1968,0.1870,0.1783,0.1706,0.1638,0.1579,0.1528,0.1484,0.1445,0.1412,0.1384,0.1359,0.1338,0.1319,0.1303,0.1288,0.1274,0.1262,0.1250,0.1239,0.1228,0.1216,0.1204,0.1191,0.1178,0.1164,0.1149,0.1133,0.1116,0.1099,0.1082,0.1064,0.1048,0.1031,0.1015,0.1000,0.0986,0.0973,0.0960,0.0948,0.0936,0.0926,0.0916,0.0907,0.0900,0.0893,0.0887,
                                                                 0.2279,0.2154,0.2038,0.1933,0.1838,0.1754,0.1680,0.1615,0.1558,0.1509,0.1466,0.1430,0.1398,0.1371,0.1348,0.1327,0.1309,0.1294,0.1280,0.1267,0.1255,0.1243,0.1232,0.1221,0.1209,0.1197,0.1184,0.1170,0.1155,0.1139,0.1123,0.1106,0.1088,0.1070,0.1052,0.1035,0.1017,0.1001,0.0985,0.0971,0.0957,0.0943,0.0931,0.0919,0.0909,0.0899,0.0890,0.0882,0.0876,0.0870,
                                                                 0.2237,0.2115,0.2003,0.1901,0.1809,0.1727,0.1656,0.1593,0.1538,0.1491,0.1451,0.1416,0.1386,0.1360,0.1337,0.1318,0.1301,0.1286,0.1272,0.1259,0.1248,0.1236,0.1225,0.1214,0.1202,0.1189,0.1176,0.1162,0.1146,0.1130,0.1113,0.1095,0.1077,0.1058,0.1039,0.1021,0.1004,0.0987,0.0970,0.0955,0.0941,0.0927,0.0915,0.0903,0.0893,0.0883,0.0875,0.0867,0.0861,0.0855,
                                                                 0.2198,0.2079,0.1970,0.1871,0.1782,0.1703,0.1634,0.1573,0.1521,0.1475,0.1436,0.1403,0.1374,0.1349,0.1328,0.1309,0.1293,0.1278,0.1265,0.1252,0.1241,0.1229,0.1218,0.1207,0.1195,0.1182,0.1168,0.1153,0.1137,0.1120,0.1102,0.1084,0.1065,0.1045,0.1026,0.1008,0.0989,0.0972,0.0956,0.0940,0.0926,0.0912,0.0900,0.0889,0.0878,0.0869,0.0861,0.0854,0.0847,0.0842,
                                                                 0.2162,0.2046,0.1939,0.1843,0.1757,0.1680,0.1613,0.1555,0.1504,0.1461,0.1423,0.1391,0.1363,0.1339,0.1319,0.1301,0.1285,0.1271,0.1258,0.1246,0.1234,0.1223,0.1211,0.1200,0.1187,0.1174,0.1160,0.1144,0.1128,0.1110,0.1091,0.1072,0.1052,0.1032,0.1013,0.0994,0.0975,0.0958,0.0941,0.0925,0.0911,0.0898,0.0886,0.0875,0.0865,0.0856,0.0848,0.0842,0.0836,0.0831,
                                                                 0.2128,0.2015,0.1911,0.1817,0.1734,0.1660,0.1595,0.1538,0.1489,0.1447,0.1411,0.1380,0.1354,0.1331,0.1311,0.1293,0.1278,0.1264,0.1251,0.1239,0.1227,0.1216,0.1204,0.1192,0.1179,0.1166,0.1151,0.1135,0.1117,0.1099,0.1080,0.1060,0.1039,0.1019,0.0999,0.0980,0.0961,0.0943,0.0927,0.0911,0.0897,0.0884,0.0873,0.0862,0.0853,0.0845,0.0838,0.0831,0.0826,0.0821,
                                                                 0.2097,0.1986,0.1885,0.1794,0.1712,0.1640,0.1577,0.1523,0.1476,0.1435,0.1400,0.1370,0.1344,0.1322,0.1303,0.1286,0.1271,0.1257,0.1244,0.1232,0.1221,0.1209,0.1197,0.1184,0.1171,0.1157,0.1141,0.1124,0.1107,0.1088,0.1068,0.1047,0.1026,0.1005,0.0985,0.0965,0.0947,0.0929,0.0913,0.0898,0.0884,0.0872,0.0861,0.0851,0.0842,0.0835,0.0828,0.0822,0.0817,0.0813,
                                                                 0.2068,0.1959,0.1861,0.1772,0.1692,0.1623,0.1562,0.1509,0.1463,0.1423,0.1390,0.1361,0.1336,0.1314,0.1295,0.1279,0.1264,0.1250,0.1238,0.1225,0.1214,0.1202,0.1189,0.1176,0.1162,0.1147,0.1131,0.1114,0.1095,0.1075,0.1055,0.1034,0.1012,0.0991,0.0971,0.0951,0.0933,0.0916,0.0900,0.0885,0.0872,0.0860,0.0850,0.0841,0.0833,0.0826,0.0820,0.0814,0.0810,0.0806,
                                                                 0.2040,0.1934,0.1838,0.1751,0.1674,0.1606,0.1547,0.1495,0.1451,0.1413,0.1380,0.1352,0.1327,0.1306,0.1288,0.1272,0.1257,0.1243,0.1231,0.1218,0.1206,0.1194,0.1181,0.1167,0.1153,0.1137,0.1120,0.1102,0.1083,0.1063,0.1042,0.1020,0.0999,0.0977,0.0957,0.0938,0.0919,0.0903,0.0887,0.0873,0.0861,0.0850,0.0840,0.0832,0.0824,0.0818,0.0812,0.0808,0.0804,0.0801,
                                                                 0.2015,0.1911,0.1817,0.1732,0.1657,0.1591,0.1533,0.1483,0.1440,0.1403,0.1371,0.1343,0.1320,0.1299,0.1281,0.1265,0.1250,0.1236,0.1223,0.1211,0.1198,0.1185,0.1172,0.1158,0.1143,0.1126,0.1109,0.1090,0.1070,0.1049,0.1028,0.1006,0.0985,0.0964,0.0943,0.0924,0.0907,0.0890,0.0876,0.0862,0.0851,0.0841,0.0832,0.0824,0.0817,0.0811,0.0806,0.0802,0.0799,0.0796,
                                                                 0.1991,0.1889,0.1797,0.1714,0.1641,0.1576,0.1520,0.1471,0.1429,0.1393,0.1362,0.1335,0.1312,0.1291,0.1274,0.1257,0.1243,0.1229,0.1216,0.1203,0.1190,0.1176,0.1162,0.1148,0.1132,0.1115,0.1097,0.1077,0.1057,0.1036,0.1014,0.0992,0.0971,0.0950,0.0930,0.0911,0.0894,0.0879,0.0865,0.0853,0.0842,0.0832,0.0824,0.0817,0.0811,0.0805,0.0801,0.0797,0.0794,0.0791,
                                                                 0.1968,0.1869,0.1778,0.1697,0.1625,0.1562,0.1508,0.1460,0.1419,0.1384,0.1353,0.1327,0.1304,0.1284,0.1266,0.1250,0.1235,0.1221,0.1208,0.1194,0.1181,0.1167,0.1152,0.1137,0.1120,0.1103,0.1084,0.1064,0.1043,0.1022,0.1000,0.0978,0.0957,0.0937,0.0917,0.0900,0.0883,0.0868,0.0855,0.0844,0.0834,0.0825,0.0817,0.0811,0.0805,0.0800,0.0796,0.0793,0.0790,0.0788,
                                                                 0.1947,0.1849,0.1761,0.1681,0.1611,0.1550,0.1496,0.1449,0.1409,0.1375,0.1345,0.1319,0.1296,0.1276,0.1259,0.1242,0.1227,0.1213,0.1199,0.1185,0.1171,0.1156,0.1141,0.1125,0.1108,0.1090,0.1071,0.1050,0.1030,0.1008,0.0986,0.0965,0.0944,0.0924,0.0906,0.0888,0.0873,0.0859,0.0846,0.0836,0.0826,0.0818,0.0811,0.0805,0.0800,0.0796,0.0793,0.0790,0.0787,0.0785,
                                                                 0.1927,0.1831,0.1744,0.1666,0.1597,0.1537,0.1485,0.1439,0.1400,0.1366,0.1336,0.1311,0.1288,0.1269,0.1251,0.1234,0.1219,0.1204,0.1190,0.1175,0.1161,0.1145,0.1130,0.1113,0.1095,0.1077,0.1057,0.1037,0.1016,0.0994,0.0973,0.0952,0.0932,0.0912,0.0895,0.0878,0.0863,0.0850,0.0839,0.0829,0.0820,0.0813,0.0806,0.0801,0.0796,0.0792,0.0789,0.0787,0.0784,0.0783,
                                                                 0.1908,0.1814,0.1728,0.1652,0.1584,0.1525,0.1474,0.1429,0.1390,0.1357,0.1328,0.1303,0.1280,0.1260,0.1242,0.1225,0.1210,0.1195,0.1180,0.1165,0.1149,0.1134,0.1117,0.1100,0.1082,0.1063,0.1044,0.1023,0.1002,0.0981,0.0960,0.0940,0.0920,0.0901,0.0884,0.0869,0.0855,0.0843,0.0832,0.0822,0.0814,0.0808,0.0802,0.0797,0.0793,0.0789,0.0787,0.0784,0.0782,0.0781,
                                                                 0.1890,0.1797,0.1713,0.1638,0.1572,0.1514,0.1463,0.1419,0.1381,0.1348,0.1319,0.1294,0.1272,0.1252,0.1233,0.1216,0.1200,0.1184,0.1169,0.1153,0.1138,0.1121,0.1105,0.1087,0.1069,0.1050,0.1030,0.1010,0.0989,0.0968,0.0948,0.0928,0.0909,0.0891,0.0875,0.0860,0.0847,0.0836,0.0826,0.0817,0.0810,0.0803,0.0798,0.0794,0.0790,0.0787,0.0784,0.0782,0.0780,0.0779,
                                                                 0.1873,0.1781,0.1698,0.1624,0.1559,0.1502,0.1452,0.1409,0.1372,0.1339,0.1310,0.1285,0.1263,0.1242,0.1224,0.1206,0.1190,0.1173,0.1157,0.1141,0.1125,0.1109,0.1091,0.1074,0.1055,0.1036,0.1017,0.0997,0.0977,0.0956,0.0937,0.0917,0.0899,0.0882,0.0867,0.0853,0.0840,0.0830,0.0820,0.0812,0.0805,0.0800,0.0795,0.0791,0.0787,0.0785,0.0782,0.0780,0.0779,0.0777,
                                                                 0.1856,0.1765,0.1684,0.1611,0.1547,0.1491,0.1442,0.1399,0.1362,0.1330,0.1301,0.1276,0.1253,0.1232,0.1213,0.1195,0.1178,0.1162,0.1145,0.1129,0.1112,0.1095,0.1078,0.1060,0.1042,0.1023,0.1004,0.0984,0.0965,0.0945,0.0926,0.0908,0.0890,0.0874,0.0859,0.0846,0.0835,0.0824,0.0816,0.0808,0.0802,0.0797,0.0792,0.0788,0.0785,0.0783,0.0781,0.0779,0.0777,0.0776,
                                                                 0.1839,0.1750,0.1670,0.1598,0.1535,0.1480,0.1431,0.1389,0.1352,0.1320,0.1291,0.1266,0.1243,0.1222,0.1202,0.1184,0.1167,0.1149,0.1133,0.1116,0.1099,0.1082,0.1065,0.1047,0.1029,0.1011,0.0992,0.0973,0.0954,0.0935,0.0916,0.0899,0.0882,0.0867,0.0853,0.0840,0.0829,0.0820,0.0812,0.0805,0.0799,0.0794,0.0790,0.0786,0.0783,0.0781,0.0779,0.0778,0.0776,0.0775,
                                                                 0.1823,0.1735,0.1656,0.1585,0.1523,0.1468,0.1420,0.1378,0.1342,0.1309,0.1281,0.1255,0.1232,0.1211,0.1191,0.1172,0.1154,0.1137,0.1120,0.1103,0.1086,0.1069,0.1052,0.1034,0.1017,0.0999,0.0980,0.0962,0.0943,0.0925,0.0908,0.0891,0.0875,0.0860,0.0847,0.0835,0.0825,0.0816,0.0808,0.0802,0.0796,0.0791,0.0788,0.0785,0.0782,0.0780,0.0778,0.0777,0.0775,0.0774,
                                                                 0.1808,0.1721,0.1642,0.1572,0.1511,0.1456,0.1409,0.1367,0.1331,0.1298,0.1270,0.1244,0.1220,0.1199,0.1178,0.1159,0.1141,0.1124,0.1106,0.1090,0.1073,0.1056,0.1039,0.1022,0.1005,0.0987,0.0970,0.0952,0.0934,0.0917,0.0900,0.0883,0.0868,0.0854,0.0842,0.0830,0.0821,0.0812,0.0805,0.0799,0.0794,0.0789,0.0786,0.0783,0.0781,0.0779,0.0777,0.0776,0.0775,0.0774,
                                                                 0.1792,0.1706,0.1628,0.1559,0.1498,0.1444,0.1397,0.1356,0.1319,0.1287,0.1258,0.1232,0.1208,0.1186,0.1166,0.1146,0.1128,0.1110,0.1093,0.1076,0.1060,0.1043,0.1027,0.1010,0.0994,0.0977,0.0960,0.0943,0.0926,0.0909,0.0893,0.0877,0.0862,0.0849,0.0837,0.0826,0.0817,0.0809,0.0802,0.0797,0.0792,0.0788,0.0784,0.0782,0.0779,0.0778,0.0776,0.0775,0.0774,0.0773,
                                                                 0.1776,0.1691,0.1614,0.1546,0.1485,0.1432,0.1385,0.1344,0.1307,0.1274,0.1245,0.1219,0.1195,0.1173,0.1152,0.1133,0.1115,0.1097,0.1080,0.1064,0.1047,0.1031,0.1016,0.1000,0.0984,0.0967,0.0951,0.0934,0.0918,0.0902,0.0886,0.0871,0.0857,0.0845,0.0833,0.0823,0.0814,0.0806,0.0800,0.0795,0.0790,0.0786,0.0783,0.0781,0.0779,0.0777,0.0775,0.0774,0.0773,0.0773,
                                                                 0.1761,0.1676,0.1600,0.1532,0.1472,0.1419,0.1372,0.1331,0.1294,0.1261,0.1232,0.1206,0.1182,0.1159,0.1139,0.1120,0.1101,0.1084,0.1067,0.1051,0.1036,0.1020,0.1005,0.0990,0.0974,0.0959,0.0943,0.0927,0.0911,0.0896,0.0881,0.0866,0.0853,0.0841,0.0829,0.0820,0.0811,0.0804,0.0798,0.0793,0.0789,0.0785,0.0782,0.0780,0.0778,0.0776,0.0775,0.0774,0.0773,0.0772,
                                                                 0.1745,0.1661,0.1585,0.1518,0.1458,0.1405,0.1359,0.1317,0.1281,0.1248,0.1219,0.1192,0.1168,0.1146,0.1125,0.1106,0.1088,0.1071,0.1055,0.1040,0.1025,0.1010,0.0995,0.0981,0.0966,0.0951,0.0936,0.0920,0.0905,0.0890,0.0876,0.0862,0.0849,0.0837,0.0826,0.0817,0.0809,0.0802,0.0796,0.0791,0.0787,0.0784,0.0781,0.0779,0.0777,0.0776,0.0774,0.0773,0.0773,0.0772,
                                                                 0.1728,0.1645,0.1570,0.1503,0.1444,0.1391,0.1344,0.1303,0.1267,0.1234,0.1205,0.1178,0.1154,0.1132,0.1112,0.1093,0.1076,0.1059,0.1044,0.1029,0.1014,0.1000,0.0986,0.0972,0.0958,0.0944,0.0929,0.0915,0.0900,0.0885,0.0871,0.0858,0.0845,0.0834,0.0824,0.0815,0.0807,0.0800,0.0795,0.0790,0.0786,0.0783,0.0780,0.0778,0.0776,0.0775,0.0774,0.0773,0.0772,0.0772,
                                                                 0.1712,0.1629,0.1554,0.1488,0.1429,0.1376,0.1330,0.1289,0.1252,0.1219,0.1190,0.1164,0.1140,0.1118,0.1099,0.1080,0.1064,0.1048,0.1033,0.1019,0.1005,0.0992,0.0978,0.0965,0.0951,0.0938,0.0924,0.0909,0.0895,0.0881,0.0867,0.0854,0.0842,0.0831,0.0821,0.0812,0.0805,0.0799,0.0793,0.0789,0.0785,0.0782,0.0780,0.0778,0.0776,0.0775,0.0774,0.0773,0.0772,0.0771,
                                                                 0.1695,0.1612,0.1538,0.1472,0.1413,0.1361,0.1314,0.1273,0.1237,0.1205,0.1176,0.1150,0.1126,0.1105,0.1086,0.1068,0.1052,0.1037,0.1023,0.1010,0.0997,0.0984,0.0971,0.0958,0.0945,0.0932,0.0918,0.0905,0.0891,0.0877,0.0864,0.0851,0.0839,0.0828,0.0819,0.0811,0.0803,0.0797,0.0792,0.0788,0.0784,0.0781,0.0779,0.0777,0.0775,0.0774,0.0773,0.0772,0.0772,0.0771,
                                                                 0.1677,0.1595,0.1521,0.1455,0.1397,0.1345,0.1299,0.1258,0.1222,0.1190,0.1161,0.1136,0.1113,0.1092,0.1074,0.1057,0.1042,0.1027,0.1014,0.1001,0.0989,0.0977,0.0965,0.0953,0.0940,0.0927,0.0914,0.0900,0.0887,0.0873,0.0860,0.0848,0.0837,0.0826,0.0817,0.0809,0.0802,0.0796,0.0791,0.0787,0.0783,0.0781,0.0778,0.0777,0.0775,0.0774,0.0773,0.0772,0.0771,0.0771,
                                                                 0.1659,0.1577,0.1504,0.1438,0.1380,0.1328,0.1283,0.1242,0.1207,0.1175,0.1147,0.1122,0.1100,0.1080,0.1063,0.1047,0.1032,0.1019,0.1006,0.0994,0.0982,0.0971,0.0959,0.0947,0.0935,0.0923,0.0910,0.0897,0.0883,0.0870,0.0857,0.0845,0.0834,0.0824,0.0815,0.0807,0.0801,0.0795,0.0790,0.0786,0.0783,0.0780,0.0778,0.0776,0.0775,0.0774,0.0773,0.0772,0.0771,0.0771,
                                                                 0.1640,0.1559,0.1486,0.1421,0.1363,0.1312,0.1266,0.1227,0.1192,0.1161,0.1134,0.1110,0.1088,0.1069,0.1052,0.1037,0.1023,0.1011,0.0999,0.0987,0.0976,0.0965,0.0954,0.0943,0.0931,0.0919,0.0906,0.0893,0.0880,0.0867,0.0855,0.0843,0.0832,0.0822,0.0814,0.0806,0.0799,0.0794,0.0789,0.0785,0.0782,0.0780,0.0778,0.0776,0.0774,0.0773,0.0772,0.0772,0.0771,0.0771};
        
        return arrayCurvatureMax[y * arrayWidth + x];
    }
    case 7:
    {
        static const std::vector<double> arrayCurvatureMax =    {0.3072,0.2950,0.2833,0.2721,0.2614,0.2512,0.2415,0.2324,0.2238,0.2158,0.2084,0.2014,0.1950,0.1891,0.1836,0.1785,0.1739,0.1696,0.1657,0.1620,0.1587,0.1555,0.1526,0.1499,0.1474,0.1449,0.1426,0.1404,0.1383,0.1362,0.1341,0.1320,0.1299,0.1277,0.1255,0.1232,0.1207,0.1182,0.1155,0.1126,0.1095,0.1062,0.1027,0.0989,0.0948,0.0904,0.0858,0.0808,0.0755,0.0699,
                                                                 0.3026,0.2906,0.2790,0.2679,0.2573,0.2472,0.2377,0.2286,0.2202,0.2123,0.2049,0.1980,0.1916,0.1858,0.1803,0.1753,0.1707,0.1665,0.1626,0.1590,0.1557,0.1527,0.1498,0.1472,0.1447,0.1424,0.1402,0.1380,0.1360,0.1340,0.1320,0.1300,0.1281,0.1261,0.1240,0.1219,0.1196,0.1173,0.1148,0.1122,0.1094,0.1064,0.1032,0.0997,0.0960,0.0920,0.0877,0.0831,0.0782,0.0730,
                                                                 0.2979,0.2859,0.2745,0.2635,0.2530,0.2431,0.2336,0.2248,0.2164,0.2086,0.2013,0.1945,0.1882,0.1823,0.1770,0.1720,0.1675,0.1633,0.1594,0.1559,0.1527,0.1497,0.1469,0.1443,0.1419,0.1396,0.1375,0.1354,0.1334,0.1315,0.1297,0.1278,0.1260,0.1241,0.1222,0.1202,0.1182,0.1161,0.1138,0.1114,0.1089,0.1061,0.1032,0.1000,0.0966,0.0930,0.0891,0.0848,0.0803,0.0755,
                                                                 0.2930,0.2812,0.2698,0.2590,0.2486,0.2388,0.2295,0.2207,0.2125,0.2047,0.1975,0.1908,0.1846,0.1788,0.1735,0.1686,0.1641,0.1599,0.1561,0.1527,0.1494,0.1465,0.1438,0.1412,0.1389,0.1367,0.1346,0.1326,0.1307,0.1289,0.1271,0.1254,0.1236,0.1219,0.1201,0.1183,0.1164,0.1145,0.1124,0.1102,0.1079,0.1054,0.1028,0.0999,0.0968,0.0935,0.0899,0.0861,0.0820,0.0775,
                                                                 0.2879,0.2762,0.2650,0.2543,0.2441,0.2344,0.2252,0.2166,0.2084,0.2008,0.1937,0.1870,0.1809,0.1752,0.1699,0.1651,0.1606,0.1565,0.1527,0.1493,0.1461,0.1432,0.1405,0.1381,0.1358,0.1336,0.1316,0.1297,0.1279,0.1261,0.1244,0.1228,0.1211,0.1195,0.1179,0.1162,0.1144,0.1126,0.1108,0.1088,0.1067,0.1044,0.1020,0.0994,0.0966,0.0936,0.0904,0.0869,0.0831,0.0790,
                                                                 0.2828,0.2712,0.2601,0.2496,0.2395,0.2299,0.2209,0.2124,0.2043,0.1968,0.1898,0.1832,0.1771,0.1715,0.1663,0.1615,0.1571,0.1530,0.1493,0.1459,0.1428,0.1399,0.1372,0.1348,0.1325,0.1304,0.1285,0.1266,0.1249,0.1232,0.1216,0.1200,0.1185,0.1169,0.1154,0.1138,0.1122,0.1106,0.1089,0.1071,0.1052,0.1031,0.1009,0.0986,0.0961,0.0934,0.0904,0.0872,0.0838,0.0801,
                                                                 0.2776,0.2661,0.2552,0.2447,0.2348,0.2254,0.2165,0.2081,0.2002,0.1928,0.1858,0.1794,0.1734,0.1678,0.1626,0.1579,0.1535,0.1495,0.1458,0.1424,0.1393,0.1365,0.1339,0.1315,0.1293,0.1272,0.1253,0.1235,0.1218,0.1202,0.1187,0.1172,0.1157,0.1143,0.1128,0.1114,0.1099,0.1084,0.1068,0.1052,0.1034,0.1016,0.0996,0.0975,0.0952,0.0928,0.0901,0.0872,0.0841,0.0808,
                                                                 0.2723,0.2610,0.2501,0.2398,0.2300,0.2207,0.2120,0.2037,0.1960,0.1887,0.1818,0.1755,0.1696,0.1641,0.1590,0.1543,0.1499,0.1460,0.1423,0.1390,0.1359,0.1331,0.1306,0.1282,0.1260,0.1240,0.1221,0.1204,0.1188,0.1172,0.1157,0.1143,0.1129,0.1115,0.1102,0.1088,0.1075,0.1061,0.1046,0.1031,0.1015,0.0999,0.0981,0.0962,0.0941,0.0919,0.0895,0.0870,0.0841,0.0811,
                                                                 0.2670,0.2558,0.2451,0.2349,0.2253,0.2161,0.2075,0.1994,0.1917,0.1846,0.1779,0.1716,0.1658,0.1604,0.1553,0.1507,0.1464,0.1425,0.1389,0.1356,0.1326,0.1298,0.1272,0.1249,0.1228,0.1208,0.1190,0.1173,0.1157,0.1142,0.1128,0.1114,0.1101,0.1088,0.1075,0.1063,0.1050,0.1037,0.1024,0.1010,0.0996,0.0980,0.0964,0.0947,0.0929,0.0909,0.0887,0.0864,0.0839,0.0812,
                                                                 0.2617,0.2507,0.2401,0.2300,0.2205,0.2115,0.2030,0.1951,0.1876,0.1805,0.1740,0.1678,0.1621,0.1567,0.1518,0.1472,0.1430,0.1391,0.1355,0.1323,0.1293,0.1265,0.1240,0.1217,0.1196,0.1177,0.1159,0.1142,0.1127,0.1112,0.1098,0.1085,0.1073,0.1061,0.1049,0.1037,0.1025,0.1013,0.1001,0.0988,0.0975,0.0962,0.0947,0.0932,0.0915,0.0897,0.0878,0.0857,0.0834,0.0810,
                                                                 0.2565,0.2456,0.2351,0.2252,0.2158,0.2070,0.1986,0.1908,0.1834,0.1765,0.1701,0.1641,0.1584,0.1532,0.1483,0.1438,0.1396,0.1358,0.1323,0.1290,0.1261,0.1233,0.1209,0.1186,0.1165,0.1146,0.1129,0.1112,0.1097,0.1083,0.1070,0.1057,0.1045,0.1034,0.1023,0.1012,0.1001,0.0990,0.0978,0.0967,0.0955,0.0942,0.0929,0.0915,0.0900,0.0884,0.0867,0.0849,0.0828,0.0806,
                                                                 0.2514,0.2405,0.2302,0.2205,0.2112,0.2025,0.1943,0.1866,0.1794,0.1726,0.1663,0.1604,0.1549,0.1497,0.1450,0.1405,0.1364,0.1326,0.1291,0.1259,0.1230,0.1203,0.1179,0.1156,0.1136,0.1117,0.1100,0.1084,0.1069,0.1055,0.1043,0.1031,0.1019,0.1008,0.0997,0.0987,0.0977,0.0966,0.0956,0.0946,0.0935,0.0923,0.0911,0.0899,0.0885,0.0871,0.0856,0.0839,0.0821,0.0801,
                                                                 0.2463,0.2356,0.2255,0.2158,0.2067,0.1982,0.1901,0.1826,0.1755,0.1689,0.1627,0.1569,0.1515,0.1464,0.1417,0.1374,0.1333,0.1296,0.1262,0.1230,0.1201,0.1174,0.1150,0.1128,0.1108,0.1089,0.1072,0.1057,0.1042,0.1029,0.1016,0.1005,0.0994,0.0983,0.0973,0.0963,0.0954,0.0944,0.0935,0.0925,0.0915,0.0905,0.0894,0.0883,0.0871,0.0858,0.0844,0.0829,0.0813,0.0796,
                                                                 0.2414,0.2309,0.2208,0.2113,0.2024,0.1940,0.1860,0.1786,0.1717,0.1652,0.1591,0.1535,0.1482,0.1433,0.1387,0.1344,0.1304,0.1268,0.1234,0.1202,0.1174,0.1147,0.1123,0.1101,0.1081,0.1063,0.1046,0.1031,0.1017,0.1004,0.0992,0.0980,0.0970,0.0960,0.0950,0.0941,0.0932,0.0923,0.0914,0.0905,0.0896,0.0887,0.0877,0.0867,0.0856,0.0845,0.0832,0.0819,0.0805,0.0789,
                                                                 0.2367,0.2263,0.2164,0.2070,0.1982,0.1899,0.1821,0.1748,0.1680,0.1617,0.1558,0.1502,0.1451,0.1402,0.1358,0.1316,0.1277,0.1241,0.1207,0.1176,0.1148,0.1122,0.1098,0.1076,0.1057,0.1038,0.1022,0.1007,0.0993,0.0980,0.0969,0.0958,0.0947,0.0938,0.0929,0.0920,0.0911,0.0903,0.0895,0.0887,0.0878,0.0870,0.0861,0.0852,0.0842,0.0832,0.0821,0.0809,0.0796,0.0782,
                                                                 0.2321,0.2218,0.2120,0.2028,0.1942,0.1860,0.1784,0.1712,0.1646,0.1583,0.1525,0.1471,0.1421,0.1374,0.1330,0.1289,0.1251,0.1216,0.1183,0.1153,0.1124,0.1099,0.1075,0.1054,0.1034,0.1016,0.1000,0.0985,0.0971,0.0959,0.0947,0.0936,0.0927,0.0917,0.0909,0.0900,0.0892,0.0884,0.0877,0.0869,0.0861,0.0854,0.0846,0.0837,0.0829,0.0819,0.0810,0.0799,0.0788,0.0776,
                                                                 0.2277,0.2175,0.2079,0.1988,0.1903,0.1823,0.1748,0.1678,0.1613,0.1552,0.1495,0.1442,0.1393,0.1347,0.1304,0.1264,0.1227,0.1192,0.1160,0.1130,0.1103,0.1077,0.1054,0.1033,0.1013,0.0995,0.0979,0.0964,0.0951,0.0939,0.0927,0.0917,0.0907,0.0898,0.0890,0.0882,0.0875,0.0867,0.0860,0.0853,0.0846,0.0839,0.0831,0.0824,0.0816,0.0808,0.0799,0.0790,0.0780,0.0769,
                                                                 0.2235,0.2135,0.2040,0.1950,0.1866,0.1788,0.1714,0.1645,0.1581,0.1522,0.1466,0.1415,0.1367,0.1322,0.1280,0.1241,0.1205,0.1171,0.1140,0.1110,0.1083,0.1058,0.1035,0.1014,0.0994,0.0977,0.0961,0.0946,0.0933,0.0921,0.0909,0.0899,0.0890,0.0881,0.0873,0.0866,0.0858,0.0851,0.0845,0.0838,0.0831,0.0825,0.0818,0.0811,0.0804,0.0797,0.0789,0.0781,0.0772,0.0762,
                                                                 0.2195,0.2096,0.2002,0.1914,0.1831,0.1754,0.1682,0.1614,0.1552,0.1493,0.1439,0.1389,0.1342,0.1298,0.1258,0.1220,0.1184,0.1151,0.1121,0.1092,0.1065,0.1041,0.1018,0.0997,0.0978,0.0960,0.0944,0.0929,0.0916,0.0904,0.0893,0.0883,0.0874,0.0866,0.0858,0.0850,0.0843,0.0837,0.0831,0.0824,0.0818,0.0812,0.0806,0.0800,0.0794,0.0787,0.0780,0.0772,0.0765,0.0756,
                                                                 0.2156,0.2059,0.1966,0.1880,0.1798,0.1722,0.1651,0.1585,0.1524,0.1467,0.1414,0.1364,0.1319,0.1276,0.1237,0.1200,0.1165,0.1133,0.1103,0.1075,0.1049,0.1025,0.1002,0.0982,0.0963,0.0945,0.0929,0.0915,0.0901,0.0890,0.0879,0.0869,0.0860,0.0851,0.0844,0.0837,0.0830,0.0824,0.0818,0.0812,0.0806,0.0801,0.0795,0.0789,0.0784,0.0778,0.0771,0.0765,0.0758,0.0750,
                                                                 0.2120,0.2023,0.1933,0.1847,0.1767,0.1692,0.1623,0.1558,0.1498,0.1442,0.1390,0.1342,0.1297,0.1256,0.1217,0.1181,0.1148,0.1117,0.1087,0.1060,0.1035,0.1011,0.0989,0.0968,0.0949,0.0932,0.0916,0.0902,0.0888,0.0877,0.0866,0.0856,0.0847,0.0839,0.0831,0.0824,0.0818,0.0812,0.0806,0.0801,0.0796,0.0790,0.0785,0.0780,0.0775,0.0769,0.0764,0.0758,0.0751,0.0745,
                                                                 0.2085,0.1990,0.1900,0.1816,0.1738,0.1664,0.1596,0.1532,0.1473,0.1418,0.1368,0.1321,0.1277,0.1237,0.1199,0.1165,0.1132,0.1102,0.1073,0.1047,0.1022,0.0998,0.0977,0.0956,0.0938,0.0920,0.0905,0.0890,0.0877,0.0865,0.0854,0.0845,0.0836,0.0828,0.0820,0.0814,0.0807,0.0801,0.0796,0.0791,0.0786,0.0781,0.0776,0.0771,0.0767,0.0762,0.0757,0.0751,0.0746,0.0740,
                                                                 0.2052,0.1958,0.1870,0.1787,0.1710,0.1637,0.1570,0.1508,0.1450,0.1396,0.1347,0.1301,0.1259,0.1219,0.1183,0.1149,0.1117,0.1088,0.1060,0.1034,0.1010,0.0987,0.0966,0.0946,0.0928,0.0910,0.0895,0.0880,0.0867,0.0855,0.0845,0.0835,0.0826,0.0818,0.0811,0.0804,0.0798,0.0792,0.0787,0.0782,0.0777,0.0773,0.0768,0.0764,0.0759,0.0755,0.0750,0.0745,0.0740,0.0735,
                                                                 0.2020,0.1928,0.1841,0.1759,0.1683,0.1612,0.1546,0.1485,0.1428,0.1376,0.1327,0.1283,0.1241,0.1203,0.1167,0.1134,0.1104,0.1075,0.1048,0.1023,0.1000,0.0977,0.0957,0.0937,0.0919,0.0902,0.0886,0.0872,0.0859,0.0847,0.0836,0.0826,0.0817,0.0809,0.0802,0.0795,0.0789,0.0784,0.0779,0.0774,0.0769,0.0765,0.0761,0.0757,0.0753,0.0749,0.0744,0.0740,0.0736,0.0731,
                                                                 0.1990,0.1899,0.1813,0.1733,0.1658,0.1588,0.1524,0.1464,0.1408,0.1357,0.1309,0.1265,0.1225,0.1188,0.1153,0.1121,0.1091,0.1064,0.1038,0.1013,0.0990,0.0969,0.0948,0.0929,0.0911,0.0895,0.0879,0.0865,0.0852,0.0840,0.0829,0.0819,0.0810,0.0802,0.0795,0.0788,0.0782,0.0777,0.0772,0.0767,0.0763,0.0759,0.0755,0.0751,0.0747,0.0743,0.0739,0.0735,0.0731,0.0727,
                                                                 0.1961,0.1871,0.1787,0.1708,0.1634,0.1566,0.1502,0.1443,0.1389,0.1339,0.1292,0.1249,0.1210,0.1174,0.1140,0.1109,0.1080,0.1053,0.1028,0.1004,0.0982,0.0961,0.0941,0.0923,0.0905,0.0889,0.0873,0.0859,0.0846,0.0834,0.0823,0.0813,0.0804,0.0796,0.0788,0.0782,0.0776,0.0770,0.0766,0.0761,0.0757,0.0753,0.0749,0.0745,0.0742,0.0738,0.0735,0.0731,0.0727,0.0724,
                                                                 0.1934,0.1845,0.1762,0.1684,0.1612,0.1545,0.1482,0.1424,0.1371,0.1322,0.1276,0.1234,0.1196,0.1160,0.1128,0.1097,0.1069,0.1043,0.1019,0.0996,0.0974,0.0954,0.0935,0.0917,0.0900,0.0883,0.0868,0.0854,0.0841,0.0829,0.0818,0.0808,0.0799,0.0791,0.0783,0.0777,0.0771,0.0765,0.0760,0.0756,0.0752,0.0748,0.0744,0.0741,0.0737,0.0734,0.0731,0.0727,0.0724,0.0721,
                                                                 0.1907,0.1820,0.1738,0.1662,0.1590,0.1524,0.1463,0.1406,0.1354,0.1306,0.1261,0.1220,0.1183,0.1148,0.1116,0.1087,0.1059,0.1034,0.1011,0.0988,0.0968,0.0948,0.0929,0.0912,0.0895,0.0879,0.0864,0.0850,0.0837,0.0825,0.0814,0.0804,0.0795,0.0786,0.0779,0.0772,0.0766,0.0761,0.0756,0.0751,0.0747,0.0743,0.0740,0.0736,0.0733,0.0730,0.0727,0.0724,0.0721,0.0718,
                                                                 0.1882,0.1796,0.1715,0.1640,0.1570,0.1505,0.1445,0.1389,0.1338,0.1291,0.1247,0.1207,0.1170,0.1137,0.1106,0.1077,0.1050,0.1026,0.1003,0.0982,0.0961,0.0942,0.0924,0.0907,0.0891,0.0876,0.0861,0.0847,0.0834,0.0822,0.0811,0.0801,0.0792,0.0783,0.0775,0.0769,0.0762,0.0757,0.0752,0.0747,0.0743,0.0740,0.0736,0.0733,0.0730,0.0727,0.0724,0.0721,0.0718,0.0716,
                                                                 0.1858,0.1773,0.1693,0.1619,0.1550,0.1487,0.1428,0.1373,0.1323,0.1276,0.1234,0.1195,0.1159,0.1126,0.1096,0.1068,0.1042,0.1018,0.0996,0.0975,0.0956,0.0937,0.0920,0.0904,0.0888,0.0873,0.0859,0.0845,0.0832,0.0820,0.0809,0.0799,0.0789,0.0781,0.0773,0.0766,0.0759,0.0754,0.0749,0.0744,0.0740,0.0736,0.0733,0.0730,0.0727,0.0724,0.0721,0.0719,0.0716,0.0714,
                                                                 0.1834,0.1751,0.1672,0.1599,0.1532,0.1469,0.1411,0.1358,0.1308,0.1263,0.1221,0.1183,0.1148,0.1116,0.1086,0.1059,0.1034,0.1011,0.0989,0.0969,0.0951,0.0933,0.0916,0.0900,0.0885,0.0871,0.0857,0.0843,0.0831,0.0819,0.0808,0.0797,0.0788,0.0779,0.0771,0.0764,0.0757,0.0751,0.0746,0.0742,0.0738,0.0734,0.0730,0.0727,0.0724,0.0722,0.0719,0.0717,0.0714,0.0712,
                                                                 0.1812,0.1729,0.1652,0.1580,0.1514,0.1452,0.1395,0.1343,0.1294,0.1250,0.1209,0.1172,0.1138,0.1106,0.1077,0.1051,0.1027,0.1004,0.0983,0.0964,0.0946,0.0929,0.0913,0.0897,0.0883,0.0869,0.0855,0.0842,0.0830,0.0818,0.0807,0.0796,0.0787,0.0778,0.0770,0.0762,0.0756,0.0750,0.0744,0.0740,0.0735,0.0732,0.0728,0.0725,0.0722,0.0720,0.0717,0.0715,0.0712,0.0710,
                                                                 0.1790,0.1709,0.1633,0.1562,0.1497,0.1436,0.1380,0.1329,0.1281,0.1238,0.1198,0.1161,0.1128,0.1097,0.1069,0.1043,0.1020,0.0998,0.0978,0.0959,0.0941,0.0925,0.0909,0.0895,0.0881,0.0867,0.0854,0.0841,0.0829,0.0818,0.0807,0.0796,0.0786,0.0777,0.0769,0.0762,0.0755,0.0749,0.0743,0.0738,0.0734,0.0730,0.0727,0.0723,0.0721,0.0718,0.0715,0.0713,0.0711,0.0709,
                                                                 0.1769,0.1689,0.1614,0.1544,0.1480,0.1420,0.1366,0.1315,0.1269,0.1226,0.1187,0.1151,0.1119,0.1089,0.1061,0.1036,0.1013,0.0992,0.0972,0.0954,0.0937,0.0921,0.0906,0.0892,0.0879,0.0866,0.0853,0.0841,0.0829,0.0818,0.0807,0.0796,0.0787,0.0778,0.0769,0.0761,0.0754,0.0748,0.0742,0.0737,0.0733,0.0729,0.0725,0.0722,0.0719,0.0717,0.0714,0.0712,0.0710,0.0708,
                                                                 0.1749,0.1670,0.1596,0.1527,0.1464,0.1405,0.1352,0.1302,0.1257,0.1215,0.1177,0.1142,0.1110,0.1081,0.1054,0.1029,0.1007,0.0986,0.0967,0.0950,0.0933,0.0918,0.0904,0.0890,0.0877,0.0864,0.0852,0.0840,0.0829,0.0818,0.0807,0.0797,0.0787,0.0778,0.0770,0.0762,0.0755,0.0748,0.0742,0.0737,0.0732,0.0728,0.0724,0.0721,0.0718,0.0715,0.0713,0.0711,0.0709,0.0707,
                                                                 0.1730,0.1651,0.1578,0.1511,0.1448,0.1391,0.1338,0.1289,0.1245,0.1204,0.1167,0.1132,0.1101,0.1073,0.1047,0.1023,0.1001,0.0981,0.0962,0.0945,0.0930,0.0915,0.0901,0.0888,0.0875,0.0863,0.0852,0.0840,0.0829,0.0819,0.0808,0.0798,0.0788,0.0779,0.0771,0.0763,0.0755,0.0748,0.0742,0.0737,0.0732,0.0728,0.0724,0.0720,0.0717,0.0715,0.0712,0.0710,0.0708,0.0706,
                                                                 0.1711,0.1633,0.1561,0.1495,0.1433,0.1377,0.1325,0.1277,0.1234,0.1194,0.1157,0.1124,0.1093,0.1065,0.1040,0.1017,0.0995,0.0976,0.0958,0.0941,0.0926,0.0912,0.0898,0.0886,0.0874,0.0862,0.0851,0.0840,0.0830,0.0819,0.0809,0.0799,0.0790,0.0781,0.0772,0.0764,0.0756,0.0749,0.0743,0.0737,0.0732,0.0728,0.0724,0.0720,0.0717,0.0714,0.0712,0.0709,0.0707,0.0705,
                                                                 0.1692,0.1616,0.1545,0.1479,0.1419,0.1363,0.1312,0.1265,0.1223,0.1183,0.1148,0.1115,0.1085,0.1058,0.1033,0.1011,0.0990,0.0971,0.0953,0.0937,0.0923,0.0909,0.0896,0.0884,0.0872,0.0861,0.0851,0.0840,0.0830,0.0820,0.0810,0.0801,0.0791,0.0782,0.0774,0.0765,0.0758,0.0751,0.0744,0.0738,0.0733,0.0728,0.0724,0.0720,0.0717,0.0714,0.0711,0.0709,0.0707,0.0705,
                                                                 0.1674,0.1599,0.1529,0.1464,0.1405,0.1350,0.1300,0.1254,0.1212,0.1174,0.1139,0.1107,0.1078,0.1051,0.1027,0.1005,0.0985,0.0966,0.0949,0.0934,0.0919,0.0906,0.0894,0.0882,0.0871,0.0860,0.0850,0.0840,0.0831,0.0821,0.0812,0.0802,0.0793,0.0784,0.0776,0.0767,0.0760,0.0752,0.0746,0.0739,0.0734,0.0729,0.0724,0.0720,0.0717,0.0714,0.0711,0.0709,0.0707,0.0705,
                                                                 0.1657,0.1583,0.1514,0.1450,0.1391,0.1337,0.1288,0.1243,0.1202,0.1164,0.1130,0.1099,0.1070,0.1044,0.1021,0.0999,0.0979,0.0962,0.0945,0.0930,0.0916,0.0903,0.0891,0.0880,0.0870,0.0859,0.0850,0.0840,0.0831,0.0822,0.0813,0.0804,0.0795,0.0786,0.0778,0.0770,0.0762,0.0754,0.0747,0.0741,0.0735,0.0730,0.0725,0.0721,0.0717,0.0714,0.0711,0.0709,0.0707,0.0705,
                                                                 0.1640,0.1567,0.1499,0.1436,0.1378,0.1325,0.1276,0.1232,0.1192,0.1155,0.1121,0.1091,0.1063,0.1038,0.1015,0.0994,0.0975,0.0957,0.0941,0.0927,0.0913,0.0901,0.0889,0.0878,0.0868,0.0859,0.0849,0.0840,0.0832,0.0823,0.0814,0.0806,0.0797,0.0789,0.0780,0.0772,0.0764,0.0757,0.0749,0.0743,0.0737,0.0731,0.0726,0.0722,0.0718,0.0715,0.0712,0.0709,0.0707,0.0705,
                                                                 0.1624,0.1551,0.1484,0.1422,0.1365,0.1313,0.1265,0.1222,0.1182,0.1146,0.1113,0.1083,0.1056,0.1031,0.1009,0.0988,0.0970,0.0953,0.0937,0.0923,0.0910,0.0898,0.0887,0.0877,0.0867,0.0858,0.0849,0.0840,0.0832,0.0824,0.0815,0.0807,0.0799,0.0791,0.0783,0.0775,0.0767,0.0759,0.0752,0.0745,0.0739,0.0733,0.0728,0.0723,0.0719,0.0716,0.0712,0.0710,0.0707,0.0705,
                                                                 0.1608,0.1536,0.1470,0.1409,0.1353,0.1301,0.1254,0.1211,0.1172,0.1137,0.1105,0.1075,0.1049,0.1025,0.1003,0.0983,0.0965,0.0949,0.0934,0.0920,0.0907,0.0896,0.0885,0.0875,0.0865,0.0857,0.0848,0.0840,0.0832,0.0824,0.0816,0.0809,0.0801,0.0793,0.0785,0.0777,0.0769,0.0762,0.0755,0.0748,0.0741,0.0735,0.0730,0.0725,0.0721,0.0717,0.0713,0.0710,0.0708,0.0705,
                                                                 0.1592,0.1521,0.1456,0.1396,0.1340,0.1290,0.1244,0.1201,0.1163,0.1128,0.1097,0.1068,0.1042,0.1019,0.0997,0.0978,0.0960,0.0944,0.0930,0.0916,0.0904,0.0893,0.0883,0.0873,0.0864,0.0856,0.0848,0.0840,0.0832,0.0825,0.0817,0.0810,0.0803,0.0795,0.0787,0.0780,0.0772,0.0765,0.0757,0.0750,0.0744,0.0738,0.0732,0.0727,0.0722,0.0718,0.0714,0.0711,0.0708,0.0706,
                                                                 0.1577,0.1507,0.1443,0.1383,0.1329,0.1279,0.1233,0.1192,0.1154,0.1120,0.1089,0.1061,0.1036,0.1013,0.0992,0.0973,0.0956,0.0940,0.0926,0.0913,0.0901,0.0891,0.0881,0.0871,0.0863,0.0855,0.0847,0.0840,0.0832,0.0825,0.0818,0.0811,0.0804,0.0797,0.0790,0.0782,0.0775,0.0768,0.0760,0.0753,0.0747,0.0740,0.0734,0.0729,0.0724,0.0720,0.0716,0.0712,0.0709,0.0707,
                                                                 0.1562,0.1493,0.1429,0.1371,0.1317,0.1268,0.1223,0.1182,0.1145,0.1112,0.1081,0.1054,0.1029,0.1007,0.0986,0.0968,0.0951,0.0936,0.0923,0.0910,0.0899,0.0888,0.0878,0.0870,0.0861,0.0853,0.0846,0.0839,0.0832,0.0826,0.0819,0.0812,0.0806,0.0799,0.0792,0.0785,0.0778,0.0771,0.0763,0.0756,0.0750,0.0743,0.0737,0.0731,0.0726,0.0722,0.0717,0.0714,0.0711,0.0708,
                                                                 0.1547,0.1479,0.1417,0.1359,0.1306,0.1258,0.1214,0.1173,0.1137,0.1104,0.1074,0.1047,0.1023,0.1001,0.0981,0.0963,0.0947,0.0932,0.0919,0.0907,0.0896,0.0886,0.0876,0.0868,0.0860,0.0852,0.0845,0.0839,0.0832,0.0826,0.0820,0.0813,0.0807,0.0801,0.0794,0.0787,0.0781,0.0774,0.0767,0.0760,0.0753,0.0746,0.0740,0.0734,0.0729,0.0724,0.0719,0.0715,0.0712,0.0709,
                                                                 0.1533,0.1466,0.1404,0.1347,0.1295,0.1248,0.1204,0.1165,0.1129,0.1096,0.1067,0.1041,0.1017,0.0995,0.0976,0.0958,0.0943,0.0928,0.0916,0.0904,0.0893,0.0883,0.0874,0.0866,0.0858,0.0851,0.0844,0.0838,0.0832,0.0826,0.0820,0.0814,0.0808,0.0802,0.0796,0.0790,0.0783,0.0777,0.0770,0.0763,0.0756,0.0749,0.0743,0.0737,0.0731,0.0726,0.0721,0.0717,0.0713,0.0710,
                                                                 0.1519,0.1453,0.1392,0.1336,0.1285,0.1238,0.1195,0.1156,0.1121,0.1089,0.1060,0.1034,0.1011,0.0990,0.0971,0.0954,0.0938,0.0925,0.0912,0.0901,0.0890,0.0881,0.0872,0.0864,0.0857,0.0850,0.0844,0.0838,0.0832,0.0826,0.0821,0.0815,0.0809,0.0804,0.0798,0.0792,0.0786,0.0779,0.0773,0.0766,0.0759,0.0753,0.0746,0.0740,0.0734,0.0729,0.0724,0.0719,0.0715,0.0712,
                                                                 0.1505,0.1440,0.1380,0.1325,0.1274,0.1228,0.1186,0.1148,0.1113,0.1082,0.1053,0.1028,0.1005,0.0984,0.0966,0.0949,0.0934,0.0921,0.0908,0.0897,0.0887,0.0878,0.0870,0.0862,0.0855,0.0849,0.0843,0.0837,0.0831,0.0826,0.0821,0.0816,0.0810,0.0805,0.0800,0.0794,0.0788,0.0782,0.0776,0.0769,0.0763,0.0756,0.0750,0.0743,0.0737,0.0732,0.0726,0.0722,0.0717,0.0714};
        
        return arrayCurvatureMax[y * arrayWidth + x];
    }
    case 8:
    {
        static const std::vector<double> arrayCurvatureMax =    {0.2979,0.2842,0.2712,0.2591,0.2478,0.2372,0.2275,0.2184,0.2101,0.2024,0.1953,0.1888,0.1828,0.1773,0.1722,0.1675,0.1632,0.1592,0.1555,0.1521,0.1490,0.1460,0.1432,0.1406,0.1382,0.1359,0.1337,0.1316,0.1296,0.1277,0.1258,0.1241,0.1223,0.1207,0.1191,0.1175,0.1159,0.1144,0.1129,0.1114,0.1100,0.1086,0.1072,0.1058,0.1044,0.1030,0.1017,0.1004,0.0990,0.0977,
                                                                 0.2924,0.2787,0.2659,0.2539,0.2428,0.2324,0.2227,0.2138,0.2056,0.1980,0.1910,0.1846,0.1787,0.1733,0.1683,0.1637,0.1595,0.1556,0.1519,0.1486,0.1454,0.1425,0.1398,0.1373,0.1348,0.1326,0.1304,0.1284,0.1264,0.1245,0.1227,0.1210,0.1193,0.1177,0.1161,0.1145,0.1130,0.1115,0.1100,0.1086,0.1072,0.1058,0.1044,0.1031,0.1017,0.1004,0.0991,0.0978,0.0965,0.0953,
                                                                 0.2868,0.2734,0.2607,0.2489,0.2378,0.2276,0.2181,0.2093,0.2012,0.1937,0.1869,0.1805,0.1747,0.1694,0.1645,0.1600,0.1558,0.1519,0.1484,0.1451,0.1420,0.1391,0.1365,0.1340,0.1316,0.1294,0.1272,0.1252,0.1233,0.1214,0.1197,0.1180,0.1163,0.1147,0.1131,0.1116,0.1101,0.1087,0.1073,0.1059,0.1045,0.1031,0.1018,0.1005,0.0992,0.0979,0.0966,0.0954,0.0941,0.0929,
                                                                 0.2814,0.2681,0.2555,0.2439,0.2330,0.2229,0.2135,0.2048,0.1968,0.1895,0.1827,0.1765,0.1708,0.1656,0.1607,0.1563,0.1522,0.1484,0.1449,0.1416,0.1386,0.1358,0.1332,0.1307,0.1284,0.1262,0.1241,0.1221,0.1203,0.1184,0.1167,0.1150,0.1134,0.1118,0.1103,0.1088,0.1074,0.1060,0.1046,0.1032,0.1019,0.1005,0.0992,0.0980,0.0967,0.0955,0.0942,0.0930,0.0918,0.0907,
                                                                 0.2760,0.2628,0.2505,0.2389,0.2282,0.2182,0.2090,0.2004,0.1926,0.1853,0.1787,0.1726,0.1670,0.1618,0.1570,0.1527,0.1486,0.1449,0.1415,0.1383,0.1353,0.1326,0.1300,0.1276,0.1253,0.1231,0.1211,0.1191,0.1173,0.1155,0.1138,0.1122,0.1106,0.1091,0.1076,0.1061,0.1047,0.1033,0.1020,0.1006,0.0993,0.0980,0.0968,0.0955,0.0943,0.0931,0.0919,0.0908,0.0896,0.0885,
                                                                 0.2707,0.2577,0.2455,0.2341,0.2235,0.2136,0.2045,0.1961,0.1884,0.1813,0.1747,0.1687,0.1632,0.1581,0.1534,0.1491,0.1452,0.1415,0.1381,0.1350,0.1321,0.1294,0.1268,0.1245,0.1222,0.1201,0.1181,0.1162,0.1144,0.1127,0.1110,0.1094,0.1079,0.1064,0.1049,0.1035,0.1021,0.1008,0.0994,0.0982,0.0969,0.0956,0.0944,0.0932,0.0920,0.0909,0.0898,0.0886,0.0875,0.0864,
                                                                 0.2655,0.2526,0.2405,0.2293,0.2188,0.2091,0.2002,0.1919,0.1843,0.1773,0.1708,0.1649,0.1595,0.1545,0.1499,0.1457,0.1418,0.1382,0.1349,0.1318,0.1290,0.1263,0.1238,0.1215,0.1193,0.1172,0.1153,0.1134,0.1116,0.1099,0.1083,0.1067,0.1052,0.1038,0.1023,0.1010,0.0996,0.0983,0.0970,0.0958,0.0945,0.0933,0.0922,0.0910,0.0899,0.0888,0.0877,0.0866,0.0855,0.0845,
                                                                 0.2603,0.2476,0.2357,0.2246,0.2143,0.2047,0.1959,0.1877,0.1803,0.1734,0.1670,0.1612,0.1559,0.1510,0.1465,0.1423,0.1385,0.1350,0.1317,0.1287,0.1259,0.1233,0.1208,0.1186,0.1164,0.1144,0.1125,0.1106,0.1089,0.1073,0.1057,0.1041,0.1027,0.1012,0.0999,0.0985,0.0972,0.0959,0.0947,0.0935,0.0923,0.0911,0.0900,0.0889,0.0878,0.0867,0.0857,0.0846,0.0836,0.0826,
                                                                 0.2553,0.2427,0.2309,0.2200,0.2098,0.2004,0.1917,0.1837,0.1763,0.1695,0.1633,0.1576,0.1523,0.1475,0.1431,0.1390,0.1353,0.1318,0.1286,0.1257,0.1229,0.1204,0.1180,0.1157,0.1136,0.1117,0.1098,0.1080,0.1063,0.1047,0.1031,0.1017,0.1002,0.0988,0.0975,0.0962,0.0949,0.0937,0.0925,0.0913,0.0902,0.0891,0.0880,0.0869,0.0858,0.0848,0.0838,0.0828,0.0819,0.0809,
                                                                 0.2503,0.2379,0.2263,0.2155,0.2055,0.1962,0.1876,0.1797,0.1725,0.1658,0.1597,0.1541,0.1489,0.1442,0.1398,0.1359,0.1322,0.1288,0.1257,0.1228,0.1201,0.1176,0.1152,0.1130,0.1110,0.1090,0.1072,0.1055,0.1038,0.1022,0.1007,0.0993,0.0979,0.0965,0.0952,0.0940,0.0927,0.0915,0.0904,0.0893,0.0881,0.0871,0.0860,0.0850,0.0840,0.0830,0.0820,0.0811,0.0802,0.0793,
                                                                 0.2454,0.2332,0.2217,0.2111,0.2012,0.1921,0.1836,0.1759,0.1687,0.1622,0.1562,0.1507,0.1456,0.1410,0.1367,0.1328,0.1292,0.1258,0.1228,0.1199,0.1173,0.1148,0.1125,0.1104,0.1084,0.1065,0.1047,0.1030,0.1014,0.0999,0.0984,0.0970,0.0956,0.0943,0.0931,0.0918,0.0907,0.0895,0.0884,0.0873,0.0862,0.0852,0.0842,0.0832,0.0822,0.0813,0.0804,0.0795,0.0786,0.0778,
                                                                 0.2407,0.2286,0.2173,0.2068,0.1971,0.1881,0.1798,0.1721,0.1651,0.1587,0.1528,0.1473,0.1424,0.1378,0.1336,0.1298,0.1263,0.1230,0.1200,0.1172,0.1146,0.1122,0.1100,0.1079,0.1059,0.1041,0.1023,0.1007,0.0991,0.0976,0.0962,0.0948,0.0935,0.0922,0.0910,0.0898,0.0887,0.0876,0.0865,0.0854,0.0844,0.0834,0.0825,0.0815,0.0806,0.0797,0.0788,0.0780,0.0772,0.0764,
                                                                 0.2360,0.2241,0.2129,0.2026,0.1930,0.1842,0.1760,0.1685,0.1616,0.1553,0.1494,0.1441,0.1393,0.1348,0.1307,0.1269,0.1234,0.1203,0.1173,0.1146,0.1121,0.1097,0.1075,0.1055,0.1036,0.1018,0.1001,0.0984,0.0969,0.0955,0.0941,0.0927,0.0915,0.0902,0.0891,0.0879,0.0868,0.0857,0.0847,0.0837,0.0827,0.0818,0.0809,0.0800,0.0791,0.0782,0.0774,0.0766,0.0758,0.0751,
                                                                 0.2315,0.2197,0.2087,0.1985,0.1891,0.1804,0.1723,0.1650,0.1582,0.1520,0.1463,0.1410,0.1362,0.1319,0.1278,0.1241,0.1208,0.1176,0.1147,0.1121,0.1096,0.1073,0.1052,0.1032,0.1013,0.0996,0.0979,0.0963,0.0948,0.0934,0.0921,0.0908,0.0896,0.0884,0.0872,0.0861,0.0851,0.0840,0.0830,0.0821,0.0811,0.0802,0.0794,0.0785,0.0777,0.0769,0.0761,0.0753,0.0746,0.0739,
                                                                 0.2270,0.2154,0.2046,0.1946,0.1853,0.1767,0.1688,0.1615,0.1549,0.1488,0.1432,0.1380,0.1334,0.1291,0.1251,0.1215,0.1182,0.1151,0.1123,0.1097,0.1073,0.1050,0.1029,0.1010,0.0992,0.0975,0.0958,0.0943,0.0929,0.0915,0.0902,0.0890,0.0878,0.0866,0.0855,0.0845,0.0834,0.0824,0.0815,0.0806,0.0797,0.0788,0.0780,0.0772,0.0764,0.0756,0.0749,0.0742,0.0735,0.0728,
                                                                 0.2227,0.2113,0.2006,0.1907,0.1816,0.1731,0.1654,0.1582,0.1517,0.1457,0.1402,0.1352,0.1306,0.1264,0.1225,0.1190,0.1157,0.1127,0.1100,0.1074,0.1050,0.1029,0.1008,0.0989,0.0971,0.0955,0.0939,0.0924,0.0910,0.0897,0.0884,0.0872,0.0861,0.0850,0.0839,0.0829,0.0819,0.0810,0.0800,0.0792,0.0783,0.0775,0.0767,0.0759,0.0752,0.0745,0.0738,0.0731,0.0724,0.0718,
                                                                 0.2186,0.2073,0.1967,0.1870,0.1780,0.1697,0.1621,0.1551,0.1486,0.1428,0.1374,0.1324,0.1279,0.1238,0.1200,0.1165,0.1134,0.1104,0.1077,0.1052,0.1029,0.1008,0.0988,0.0970,0.0952,0.0936,0.0921,0.0907,0.0893,0.0880,0.0868,0.0856,0.0845,0.0834,0.0824,0.0814,0.0805,0.0796,0.0787,0.0779,0.0771,0.0763,0.0755,0.0748,0.0741,0.0734,0.0728,0.0721,0.0715,0.0709,
                                                                 0.2145,0.2034,0.1930,0.1834,0.1746,0.1664,0.1589,0.1520,0.1457,0.1399,0.1347,0.1298,0.1254,0.1214,0.1176,0.1143,0.1111,0.1083,0.1056,0.1032,0.1010,0.0989,0.0969,0.0951,0.0935,0.0919,0.0904,0.0890,0.0877,0.0864,0.0853,0.0841,0.0831,0.0820,0.0810,0.0801,0.0792,0.0783,0.0775,0.0767,0.0759,0.0752,0.0745,0.0738,0.0731,0.0725,0.0718,0.0713,0.0707,0.0702,
                                                                 0.2106,0.1996,0.1894,0.1800,0.1712,0.1632,0.1559,0.1491,0.1429,0.1372,0.1321,0.1273,0.1230,0.1190,0.1154,0.1121,0.1090,0.1062,0.1037,0.1013,0.0991,0.0971,0.0952,0.0934,0.0918,0.0903,0.0888,0.0875,0.0862,0.0850,0.0838,0.0827,0.0817,0.0807,0.0798,0.0789,0.0780,0.0772,0.0764,0.0756,0.0749,0.0742,0.0735,0.0728,0.0722,0.0716,0.0710,0.0705,0.0700,0.0695,
                                                                 0.2067,0.1959,0.1859,0.1766,0.1681,0.1602,0.1530,0.1463,0.1402,0.1347,0.1296,0.1250,0.1207,0.1168,0.1133,0.1100,0.1071,0.1043,0.1018,0.0995,0.0974,0.0954,0.0935,0.0918,0.0902,0.0887,0.0873,0.0860,0.0848,0.0836,0.0825,0.0815,0.0805,0.0795,0.0786,0.0778,0.0769,0.0761,0.0754,0.0746,0.0739,0.0733,0.0726,0.0720,0.0714,0.0709,0.0703,0.0698,0.0693,0.0689,
                                                                 0.2031,0.1924,0.1826,0.1734,0.1650,0.1573,0.1502,0.1437,0.1377,0.1323,0.1273,0.1227,0.1186,0.1148,0.1113,0.1081,0.1052,0.1025,0.1001,0.0978,0.0957,0.0938,0.0920,0.0903,0.0888,0.0873,0.0860,0.0847,0.0835,0.0824,0.0813,0.0803,0.0794,0.0784,0.0776,0.0767,0.0759,0.0752,0.0745,0.0738,0.0731,0.0725,0.0719,0.0713,0.0707,0.0702,0.0697,0.0692,0.0688,0.0683,
                                                                 0.1995,0.1891,0.1793,0.1704,0.1621,0.1545,0.1475,0.1411,0.1353,0.1299,0.1251,0.1206,0.1165,0.1128,0.1094,0.1063,0.1035,0.1009,0.0985,0.0963,0.0942,0.0923,0.0906,0.0890,0.0875,0.0861,0.0848,0.0835,0.0824,0.0813,0.0802,0.0793,0.0783,0.0775,0.0766,0.0758,0.0751,0.0743,0.0736,0.0730,0.0723,0.0717,0.0712,0.0706,0.0701,0.0696,0.0691,0.0687,0.0683,0.0679,
                                                                 0.1961,0.1858,0.1763,0.1675,0.1593,0.1519,0.1450,0.1387,0.1330,0.1278,0.1230,0.1186,0.1146,0.1110,0.1077,0.1046,0.1019,0.0993,0.0970,0.0948,0.0928,0.0910,0.0893,0.0877,0.0863,0.0849,0.0836,0.0824,0.0813,0.0802,0.0793,0.0783,0.0774,0.0766,0.0758,0.0750,0.0743,0.0736,0.0729,0.0723,0.0717,0.0711,0.0706,0.0700,0.0696,0.0691,0.0686,0.0682,0.0678,0.0675,
                                                                 0.1929,0.1827,0.1733,0.1647,0.1567,0.1494,0.1426,0.1365,0.1309,0.1257,0.1211,0.1168,0.1129,0.1093,0.1061,0.1031,0.1004,0.0979,0.0956,0.0935,0.0916,0.0898,0.0881,0.0866,0.0852,0.0838,0.0826,0.0814,0.0804,0.0793,0.0784,0.0775,0.0766,0.0758,0.0750,0.0743,0.0736,0.0729,0.0723,0.0717,0.0711,0.0706,0.0700,0.0695,0.0691,0.0686,0.0682,0.0678,0.0675,0.0671,
                                                                 0.1897,0.1798,0.1705,0.1620,0.1542,0.1470,0.1404,0.1344,0.1289,0.1238,0.1192,0.1151,0.1112,0.1078,0.1046,0.1017,0.0990,0.0966,0.0943,0.0923,0.0904,0.0887,0.0871,0.0856,0.0842,0.0829,0.0817,0.0806,0.0795,0.0785,0.0776,0.0767,0.0759,0.0751,0.0743,0.0736,0.0729,0.0723,0.0717,0.0711,0.0706,0.0701,0.0696,0.0691,0.0687,0.0683,0.0679,0.0675,0.0672,0.0668,
                                                                 0.1867,0.1769,0.1679,0.1595,0.1518,0.1447,0.1383,0.1324,0.1270,0.1220,0.1175,0.1135,0.1097,0.1063,0.1032,0.1004,0.0978,0.0954,0.0932,0.0912,0.0894,0.0877,0.0861,0.0846,0.0833,0.0820,0.0809,0.0798,0.0787,0.0778,0.0769,0.0760,0.0752,0.0745,0.0737,0.0731,0.0724,0.0718,0.0712,0.0707,0.0701,0.0696,0.0692,0.0687,0.0683,0.0679,0.0676,0.0672,0.0669,0.0666,
                                                                 0.1839,0.1742,0.1653,0.1571,0.1496,0.1426,0.1363,0.1305,0.1252,0.1204,0.1160,0.1120,0.1083,0.1050,0.1019,0.0992,0.0966,0.0943,0.0922,0.0902,0.0884,0.0868,0.0852,0.0838,0.0825,0.0813,0.0801,0.0791,0.0781,0.0771,0.0763,0.0754,0.0746,0.0739,0.0732,0.0726,0.0719,0.0713,0.0708,0.0702,0.0697,0.0693,0.0688,0.0684,0.0680,0.0676,0.0673,0.0670,0.0667,0.0664,
                                                                 0.1811,0.1717,0.1629,0.1548,0.1474,0.1407,0.1344,0.1288,0.1236,0.1188,0.1145,0.1106,0.1070,0.1038,0.1008,0.0981,0.0956,0.0933,0.0912,0.0893,0.0876,0.0859,0.0845,0.0831,0.0818,0.0806,0.0795,0.0784,0.0775,0.0766,0.0757,0.0749,0.0741,0.0734,0.0728,0.0721,0.0715,0.0709,0.0704,0.0699,0.0694,0.0689,0.0685,0.0681,0.0677,0.0674,0.0670,0.0667,0.0664,0.0662,
                                                                 0.1785,0.1692,0.1606,0.1527,0.1455,0.1388,0.1327,0.1271,0.1221,0.1174,0.1132,0.1094,0.1059,0.1027,0.0998,0.0971,0.0947,0.0924,0.0904,0.0885,0.0868,0.0852,0.0838,0.0824,0.0812,0.0800,0.0789,0.0779,0.0770,0.0761,0.0752,0.0744,0.0737,0.0730,0.0724,0.0717,0.0711,0.0706,0.0701,0.0696,0.0691,0.0687,0.0682,0.0679,0.0675,0.0671,0.0668,0.0665,0.0662,0.0660,
                                                                 0.1760,0.1669,0.1585,0.1507,0.1436,0.1371,0.1311,0.1256,0.1207,0.1161,0.1120,0.1082,0.1048,0.1017,0.0988,0.0962,0.0938,0.0917,0.0897,0.0878,0.0862,0.0846,0.0832,0.0818,0.0806,0.0795,0.0784,0.0774,0.0765,0.0756,0.0748,0.0741,0.0733,0.0727,0.0720,0.0714,0.0708,0.0703,0.0698,0.0693,0.0688,0.0684,0.0680,0.0676,0.0673,0.0669,0.0666,0.0663,0.0661,0.0658,
                                                                 0.1736,0.1647,0.1564,0.1488,0.1418,0.1354,0.1296,0.1242,0.1194,0.1149,0.1109,0.1072,0.1038,0.1008,0.0980,0.0954,0.0931,0.0910,0.0890,0.0872,0.0856,0.0840,0.0826,0.0813,0.0801,0.0790,0.0780,0.0770,0.0761,0.0753,0.0745,0.0737,0.0730,0.0723,0.0717,0.0711,0.0705,0.0700,0.0695,0.0690,0.0686,0.0682,0.0678,0.0674,0.0670,0.0667,0.0664,0.0661,0.0659,0.0656,
                                                                 0.1714,0.1626,0.1545,0.1470,0.1402,0.1339,0.1282,0.1229,0.1182,0.1138,0.1099,0.1062,0.1029,0.0999,0.0972,0.0947,0.0924,0.0903,0.0884,0.0867,0.0851,0.0836,0.0822,0.0809,0.0797,0.0786,0.0776,0.0767,0.0758,0.0749,0.0741,0.0734,0.0727,0.0720,0.0714,0.0708,0.0703,0.0698,0.0693,0.0688,0.0684,0.0679,0.0675,0.0672,0.0668,0.0665,0.0662,0.0659,0.0657,0.0654,
                                                                 0.1692,0.1606,0.1526,0.1453,0.1386,0.1325,0.1269,0.1218,0.1171,0.1128,0.1089,0.1054,0.1022,0.0992,0.0965,0.0941,0.0918,0.0898,0.0879,0.0862,0.0846,0.0831,0.0818,0.0805,0.0794,0.0783,0.0773,0.0763,0.0755,0.0746,0.0739,0.0731,0.0724,0.0718,0.0712,0.0706,0.0700,0.0695,0.0690,0.0686,0.0681,0.0677,0.0673,0.0669,0.0666,0.0663,0.0660,0.0657,0.0654,0.0652,
                                                                 0.1671,0.1587,0.1509,0.1438,0.1372,0.1312,0.1257,0.1206,0.1161,0.1119,0.1081,0.1046,0.1014,0.0986,0.0959,0.0935,0.0913,0.0893,0.0875,0.0858,0.0842,0.0828,0.0814,0.0802,0.0791,0.0780,0.0770,0.0761,0.0752,0.0744,0.0736,0.0729,0.0722,0.0715,0.0709,0.0704,0.0698,0.0693,0.0688,0.0683,0.0679,0.0675,0.0671,0.0667,0.0664,0.0660,0.0657,0.0655,0.0652,0.0650,
                                                                 0.1651,0.1569,0.1493,0.1423,0.1358,0.1299,0.1245,0.1196,0.1151,0.1110,0.1073,0.1039,0.1008,0.0980,0.0954,0.0930,0.0909,0.0889,0.0871,0.0854,0.0839,0.0824,0.0811,0.0799,0.0788,0.0777,0.0767,0.0758,0.0749,0.0741,0.0734,0.0726,0.0719,0.0713,0.0707,0.0701,0.0696,0.0690,0.0686,0.0681,0.0676,0.0672,0.0668,0.0665,0.0661,0.0658,0.0655,0.0652,0.0649,0.0647,
                                                                 0.1632,0.1551,0.1477,0.1408,0.1345,0.1288,0.1235,0.1187,0.1143,0.1103,0.1066,0.1033,0.1002,0.0974,0.0949,0.0926,0.0905,0.0885,0.0867,0.0851,0.0835,0.0821,0.0808,0.0796,0.0785,0.0775,0.0765,0.0756,0.0747,0.0739,0.0731,0.0724,0.0717,0.0711,0.0705,0.0699,0.0693,0.0688,0.0683,0.0678,0.0674,0.0670,0.0665,0.0662,0.0658,0.0655,0.0652,0.0649,0.0646,0.0643,
                                                                 0.1614,0.1535,0.1462,0.1395,0.1333,0.1277,0.1225,0.1178,0.1135,0.1095,0.1060,0.1027,0.0997,0.0970,0.0945,0.0922,0.0901,0.0882,0.0864,0.0848,0.0833,0.0819,0.0806,0.0794,0.0783,0.0772,0.0762,0.0753,0.0745,0.0736,0.0729,0.0722,0.0715,0.0708,0.0702,0.0696,0.0691,0.0685,0.0680,0.0675,0.0671,0.0666,0.0662,0.0658,0.0655,0.0651,0.0648,0.0645,0.0642,0.0639,
                                                                 0.1596,0.1519,0.1447,0.1382,0.1321,0.1266,0.1216,0.1169,0.1127,0.1089,0.1054,0.1021,0.0992,0.0965,0.0941,0.0918,0.0898,0.0879,0.0861,0.0845,0.0830,0.0816,0.0803,0.0791,0.0780,0.0770,0.0760,0.0751,0.0742,0.0734,0.0726,0.0719,0.0712,0.0706,0.0699,0.0693,0.0688,0.0682,0.0677,0.0672,0.0667,0.0663,0.0659,0.0655,0.0651,0.0647,0.0644,0.0641,0.0638,0.0635,
                                                                 0.1578,0.1503,0.1433,0.1369,0.1310,0.1256,0.1207,0.1162,0.1120,0.1082,0.1048,0.1016,0.0987,0.0961,0.0937,0.0915,0.0894,0.0876,0.0858,0.0842,0.0827,0.0814,0.0801,0.0789,0.0778,0.0767,0.0758,0.0748,0.0740,0.0731,0.0724,0.0716,0.0709,0.0703,0.0696,0.0690,0.0684,0.0679,0.0673,0.0668,0.0664,0.0659,0.0655,0.0651,0.0647,0.0643,0.0639,0.0636,0.0633,0.0630,
                                                                 0.1561,0.1487,0.1419,0.1357,0.1299,0.1246,0.1198,0.1154,0.1113,0.1076,0.1042,0.1011,0.0983,0.0957,0.0933,0.0911,0.0891,0.0873,0.0855,0.0840,0.0825,0.0811,0.0798,0.0787,0.0775,0.0765,0.0755,0.0746,0.0737,0.0729,0.0721,0.0713,0.0706,0.0699,0.0693,0.0687,0.0681,0.0675,0.0670,0.0664,0.0659,0.0655,0.0650,0.0646,0.0642,0.0638,0.0634,0.0631,0.0627,0.0624,
                                                                 0.1544,0.1472,0.1406,0.1345,0.1289,0.1237,0.1190,0.1146,0.1107,0.1070,0.1037,0.1007,0.0979,0.0953,0.0930,0.0908,0.0888,0.0870,0.0853,0.0837,0.0822,0.0809,0.0796,0.0784,0.0773,0.0762,0.0752,0.0743,0.0734,0.0725,0.0717,0.0710,0.0703,0.0696,0.0689,0.0683,0.0676,0.0671,0.0665,0.0660,0.0654,0.0650,0.0645,0.0640,0.0636,0.0632,0.0628,0.0624,0.0621,0.0618,
                                                                 0.1527,0.1457,0.1392,0.1333,0.1278,0.1228,0.1181,0.1139,0.1100,0.1064,0.1032,0.1002,0.0974,0.0949,0.0926,0.0905,0.0885,0.0867,0.0850,0.0834,0.0819,0.0806,0.0793,0.0781,0.0770,0.0759,0.0749,0.0740,0.0730,0.0722,0.0714,0.0706,0.0698,0.0691,0.0685,0.0678,0.0672,0.0666,0.0660,0.0654,0.0649,0.0644,0.0639,0.0634,0.0630,0.0626,0.0621,0.0618,0.0614,0.0610,
                                                                 0.1510,0.1442,0.1379,0.1321,0.1267,0.1218,0.1173,0.1131,0.1093,0.1059,0.1026,0.0997,0.0970,0.0945,0.0922,0.0901,0.0881,0.0863,0.0847,0.0831,0.0816,0.0803,0.0790,0.0778,0.0766,0.0756,0.0745,0.0736,0.0727,0.0718,0.0710,0.0702,0.0694,0.0687,0.0680,0.0673,0.0666,0.0660,0.0654,0.0648,0.0643,0.0638,0.0632,0.0628,0.0623,0.0618,0.0614,0.0610,0.0606,0.0602,
                                                                 0.1492,0.1426,0.1365,0.1308,0.1256,0.1208,0.1164,0.1124,0.1087,0.1052,0.1021,0.0992,0.0965,0.0941,0.0918,0.0897,0.0878,0.0860,0.0843,0.0827,0.0813,0.0799,0.0786,0.0774,0.0762,0.0752,0.0741,0.0732,0.0722,0.0713,0.0705,0.0697,0.0689,0.0681,0.0674,0.0667,0.0660,0.0654,0.0648,0.0642,0.0636,0.0630,0.0625,0.0620,0.0615,0.0610,0.0606,0.0602,0.0597,0.0593,
                                                                 0.1474,0.1410,0.1351,0.1296,0.1245,0.1198,0.1155,0.1116,0.1079,0.1046,0.1015,0.0987,0.0960,0.0936,0.0914,0.0893,0.0874,0.0856,0.0839,0.0823,0.0809,0.0795,0.0782,0.0770,0.0758,0.0747,0.0737,0.0727,0.0717,0.0708,0.0699,0.0691,0.0683,0.0675,0.0668,0.0661,0.0654,0.0647,0.0641,0.0634,0.0628,0.0623,0.0617,0.0612,0.0607,0.0602,0.0597,0.0592,0.0588,0.0584,
                                                                 0.1456,0.1394,0.1336,0.1283,0.1233,0.1188,0.1146,0.1107,0.1072,0.1039,0.1009,0.0981,0.0955,0.0931,0.0909,0.0888,0.0869,0.0851,0.0835,0.0819,0.0804,0.0790,0.0777,0.0765,0.0753,0.0742,0.0731,0.0721,0.0712,0.0702,0.0693,0.0685,0.0676,0.0668,0.0661,0.0653,0.0646,0.0639,0.0633,0.0626,0.0620,0.0614,0.0608,0.0603,0.0597,0.0592,0.0587,0.0582,0.0578,0.0573,
                                                                 0.1436,0.1376,0.1320,0.1269,0.1221,0.1177,0.1136,0.1099,0.1064,0.1032,0.1002,0.0975,0.0949,0.0926,0.0904,0.0883,0.0864,0.0846,0.0830,0.0814,0.0799,0.0785,0.0772,0.0760,0.0748,0.0736,0.0726,0.0715,0.0705,0.0696,0.0687,0.0678,0.0669,0.0661,0.0653,0.0645,0.0638,0.0631,0.0624,0.0617,0.0611,0.0605,0.0599,0.0593,0.0587,0.0582,0.0576,0.0571,0.0567,0.0562,
                                                                 0.1416,0.1358,0.1304,0.1254,0.1208,0.1165,0.1126,0.1089,0.1055,0.1024,0.0995,0.0968,0.0943,0.0919,0.0898,0.0877,0.0858,0.0841,0.0824,0.0808,0.0793,0.0779,0.0766,0.0753,0.0741,0.0730,0.0719,0.0708,0.0698,0.0688,0.0679,0.0670,0.0661,0.0653,0.0645,0.0637,0.0629,0.0622,0.0614,0.0608,0.0601,0.0594,0.0588,0.0582,0.0576,0.0570,0.0565,0.0560,0.0555,0.0550,
                                                                 0.1395,0.1339,0.1287,0.1239,0.1194,0.1153,0.1114,0.1079,0.1046,0.1015,0.0986,0.0960,0.0935,0.0912,0.0891,0.0871,0.0852,0.0834,0.0818,0.0802,0.0787,0.0773,0.0759,0.0747,0.0734,0.0723,0.0712,0.0701,0.0690,0.0680,0.0671,0.0661,0.0652,0.0644,0.0635,0.0627,0.0619,0.0612,0.0604,0.0597,0.0590,0.0583,0.0577,0.0570,0.0564,0.0558,0.0553,0.0547,0.0542,0.0537,
                                                                 0.1372,0.1319,0.1269,0.1222,0.1179,0.1139,0.1102,0.1067,0.1035,0.1005,0.0978,0.0952,0.0927,0.0905,0.0884,0.0864,0.0845,0.0827,0.0811,0.0795,0.0780,0.0766,0.0752,0.0739,0.0727,0.0715,0.0703,0.0692,0.0682,0.0672,0.0662,0.0652,0.0643,0.0634,0.0625,0.0617,0.0609,0.0601,0.0593,0.0586,0.0578,0.0571,0.0565,0.0558,0.0552,0.0546,0.0540,0.0534,0.0528,0.0523};
        
        return arrayCurvatureMax[y * arrayWidth + x];
    }
    case 9:
    {
        static const std::vector<double> arrayCurvatureMax =    {0.2531,0.2445,0.2359,0.2275,0.2192,0.2113,0.2037,0.1967,0.1901,0.1840,0.1784,0.1733,0.1687,0.1644,0.1605,0.1570,0.1537,0.1507,0.1479,0.1452,0.1427,0.1403,0.1379,0.1356,0.1334,0.1312,0.1290,0.1269,0.1247,0.1226,0.1204,0.1183,0.1161,0.1140,0.1119,0.1098,0.1078,0.1058,0.1038,0.1020,0.1003,0.0986,0.0972,0.0959,0.0948,0.0938,0.0931,0.0925,0.0920,0.0916,
                                                                 0.2514,0.2428,0.2342,0.2256,0.2173,0.2094,0.2018,0.1947,0.1881,0.1820,0.1764,0.1712,0.1665,0.1622,0.1582,0.1546,0.1513,0.1481,0.1452,0.1425,0.1398,0.1373,0.1349,0.1325,0.1302,0.1279,0.1257,0.1235,0.1212,0.1191,0.1169,0.1147,0.1126,0.1104,0.1083,0.1063,0.1042,0.1023,0.1004,0.0986,0.0969,0.0954,0.0940,0.0928,0.0918,0.0909,0.0902,0.0897,0.0894,0.0891,
                                                                 0.2496,0.2409,0.2322,0.2237,0.2153,0.2073,0.1997,0.1926,0.1859,0.1798,0.1741,0.1689,0.1641,0.1597,0.1557,0.1520,0.1486,0.1454,0.1424,0.1395,0.1368,0.1342,0.1317,0.1293,0.1269,0.1245,0.1222,0.1200,0.1178,0.1155,0.1134,0.1112,0.1091,0.1070,0.1049,0.1029,0.1009,0.0990,0.0972,0.0955,0.0939,0.0924,0.0911,0.0900,0.0890,0.0883,0.0877,0.0873,0.0870,0.0868,
                                                                 0.2477,0.2389,0.2302,0.2215,0.2132,0.2051,0.1975,0.1903,0.1836,0.1774,0.1716,0.1664,0.1615,0.1571,0.1530,0.1492,0.1457,0.1424,0.1393,0.1364,0.1336,0.1309,0.1284,0.1259,0.1235,0.1211,0.1188,0.1165,0.1143,0.1121,0.1099,0.1078,0.1057,0.1036,0.1016,0.0996,0.0977,0.0959,0.0942,0.0925,0.0910,0.0896,0.0884,0.0873,0.0865,0.0858,0.0853,0.0850,0.0848,0.0847,
                                                                 0.2457,0.2368,0.2279,0.2192,0.2108,0.2027,0.1950,0.1877,0.1810,0.1747,0.1689,0.1636,0.1587,0.1542,0.1500,0.1462,0.1426,0.1392,0.1361,0.1331,0.1303,0.1276,0.1249,0.1224,0.1200,0.1176,0.1153,0.1131,0.1109,0.1087,0.1066,0.1045,0.1025,0.1005,0.0985,0.0966,0.0948,0.0930,0.0914,0.0898,0.0884,0.0871,0.0859,0.0850,0.0842,0.0836,0.0832,0.0829,0.0828,0.0828,
                                                                 0.2434,0.2345,0.2255,0.2168,0.2082,0.2001,0.1923,0.1850,0.1782,0.1719,0.1660,0.1606,0.1557,0.1511,0.1469,0.1430,0.1393,0.1359,0.1327,0.1297,0.1268,0.1241,0.1215,0.1190,0.1166,0.1142,0.1119,0.1097,0.1076,0.1055,0.1034,0.1014,0.0994,0.0975,0.0956,0.0938,0.0920,0.0904,0.0888,0.0873,0.0859,0.0847,0.0837,0.0828,0.0821,0.0816,0.0812,0.0810,0.0809,0.0810,
                                                                 0.2410,0.2319,0.2229,0.2141,0.2055,0.1972,0.1894,0.1821,0.1752,0.1688,0.1629,0.1575,0.1525,0.1478,0.1436,0.1396,0.1359,0.1325,0.1293,0.1262,0.1234,0.1207,0.1181,0.1156,0.1132,0.1109,0.1087,0.1065,0.1044,0.1024,0.1004,0.0985,0.0966,0.0947,0.0929,0.0912,0.0895,0.0879,0.0864,0.0850,0.0838,0.0826,0.0816,0.0808,0.0802,0.0797,0.0794,0.0792,0.0792,0.0793,
                                                                 0.2384,0.2292,0.2201,0.2111,0.2025,0.1942,0.1863,0.1789,0.1720,0.1655,0.1596,0.1541,0.1491,0.1444,0.1401,0.1361,0.1324,0.1290,0.1258,0.1228,0.1199,0.1172,0.1147,0.1123,0.1099,0.1077,0.1056,0.1035,0.1015,0.0995,0.0976,0.0958,0.0939,0.0922,0.0905,0.0888,0.0872,0.0857,0.0843,0.0830,0.0818,0.0807,0.0798,0.0790,0.0785,0.0780,0.0778,0.0776,0.0777,0.0778,
                                                                 0.2355,0.2263,0.2170,0.2080,0.1993,0.1909,0.1830,0.1755,0.1686,0.1621,0.1561,0.1506,0.1455,0.1409,0.1366,0.1326,0.1289,0.1255,0.1223,0.1193,0.1165,0.1139,0.1114,0.1091,0.1068,0.1047,0.1026,0.1006,0.0987,0.0968,0.0950,0.0932,0.0915,0.0898,0.0882,0.0866,0.0851,0.0837,0.0823,0.0811,0.0800,0.0790,0.0781,0.0774,0.0769,0.0765,0.0763,0.0762,0.0762,0.0763,
                                                                 0.2325,0.2231,0.2138,0.2047,0.1958,0.1874,0.1794,0.1719,0.1650,0.1585,0.1525,0.1470,0.1419,0.1373,0.1330,0.1290,0.1254,0.1220,0.1189,0.1160,0.1133,0.1107,0.1083,0.1060,0.1039,0.1018,0.0999,0.0980,0.0961,0.0944,0.0926,0.0909,0.0893,0.0877,0.0861,0.0846,0.0832,0.0818,0.0806,0.0794,0.0783,0.0774,0.0766,0.0760,0.0755,0.0751,0.0749,0.0748,0.0748,0.0750,
                                                                 0.2292,0.2197,0.2103,0.2011,0.1922,0.1838,0.1757,0.1682,0.1612,0.1548,0.1488,0.1433,0.1383,0.1336,0.1294,0.1255,0.1219,0.1186,0.1156,0.1128,0.1101,0.1077,0.1054,0.1032,0.1012,0.0992,0.0973,0.0955,0.0938,0.0921,0.0904,0.0888,0.0873,0.0857,0.0843,0.0828,0.0815,0.0802,0.0790,0.0779,0.0769,0.0760,0.0752,0.0746,0.0742,0.0738,0.0736,0.0735,0.0736,0.0737,
                                                                 0.2257,0.2161,0.2066,0.1973,0.1884,0.1799,0.1719,0.1644,0.1574,0.1510,0.1450,0.1396,0.1346,0.1300,0.1259,0.1221,0.1186,0.1154,0.1124,0.1097,0.1072,0.1048,0.1026,0.1006,0.0986,0.0968,0.0950,0.0933,0.0916,0.0900,0.0885,0.0869,0.0854,0.0840,0.0825,0.0812,0.0799,0.0786,0.0775,0.0765,0.0755,0.0747,0.0740,0.0734,0.0730,0.0726,0.0724,0.0724,0.0724,0.0725,
                                                                 0.2220,0.2123,0.2027,0.1934,0.1845,0.1760,0.1680,0.1605,0.1535,0.1471,0.1413,0.1359,0.1310,0.1265,0.1225,0.1188,0.1154,0.1123,0.1095,0.1068,0.1044,0.1022,0.1001,0.0982,0.0963,0.0946,0.0929,0.0913,0.0897,0.0881,0.0866,0.0852,0.0837,0.0823,0.0810,0.0797,0.0784,0.0773,0.0762,0.0752,0.0743,0.0735,0.0729,0.0723,0.0719,0.0715,0.0713,0.0713,0.0713,0.0714,
                                                                 0.2181,0.2083,0.1987,0.1894,0.1804,0.1720,0.1640,0.1565,0.1497,0.1433,0.1376,0.1323,0.1275,0.1231,0.1192,0.1156,0.1124,0.1094,0.1067,0.1042,0.1019,0.0998,0.0978,0.0960,0.0942,0.0925,0.0909,0.0894,0.0879,0.0864,0.0850,0.0836,0.0822,0.0809,0.0796,0.0783,0.0771,0.0760,0.0750,0.0740,0.0732,0.0725,0.0718,0.0713,0.0708,0.0705,0.0703,0.0702,0.0702,0.0703,
                                                                 0.2140,0.2042,0.1946,0.1853,0.1763,0.1679,0.1600,0.1526,0.1459,0.1396,0.1340,0.1288,0.1242,0.1199,0.1161,0.1127,0.1095,0.1067,0.1041,0.1018,0.0996,0.0976,0.0957,0.0940,0.0923,0.0907,0.0892,0.0877,0.0863,0.0849,0.0835,0.0822,0.0808,0.0795,0.0783,0.0771,0.0760,0.0749,0.0739,0.0730,0.0722,0.0715,0.0709,0.0703,0.0699,0.0696,0.0694,0.0692,0.0692,0.0693,
                                                                 0.2098,0.2000,0.1904,0.1811,0.1722,0.1639,0.1561,0.1488,0.1421,0.1361,0.1305,0.1255,0.1210,0.1169,0.1132,0.1099,0.1069,0.1042,0.1018,0.0995,0.0975,0.0956,0.0938,0.0922,0.0906,0.0891,0.0876,0.0862,0.0848,0.0835,0.0822,0.0808,0.0796,0.0783,0.0771,0.0760,0.0749,0.0739,0.0729,0.0720,0.0713,0.0706,0.0700,0.0694,0.0690,0.0687,0.0684,0.0683,0.0683,0.0683,
                                                                 0.2056,0.1958,0.1862,0.1770,0.1682,0.1599,0.1522,0.1451,0.1386,0.1326,0.1273,0.1224,0.1180,0.1141,0.1106,0.1074,0.1046,0.1020,0.0997,0.0975,0.0956,0.0938,0.0921,0.0905,0.0890,0.0876,0.0862,0.0848,0.0835,0.0822,0.0809,0.0796,0.0784,0.0772,0.0760,0.0749,0.0739,0.0729,0.0720,0.0712,0.0704,0.0697,0.0691,0.0686,0.0682,0.0678,0.0676,0.0674,0.0673,0.0674,
                                                                 0.2014,0.1916,0.1820,0.1729,0.1642,0.1561,0.1485,0.1416,0.1352,0.1294,0.1242,0.1195,0.1153,0.1115,0.1081,0.1051,0.1024,0.0999,0.0977,0.0957,0.0939,0.0922,0.0906,0.0891,0.0876,0.0862,0.0849,0.0836,0.0823,0.0810,0.0798,0.0785,0.0773,0.0762,0.0750,0.0740,0.0730,0.0720,0.0712,0.0704,0.0696,0.0690,0.0684,0.0679,0.0674,0.0670,0.0668,0.0666,0.0665,0.0665,
                                                                 0.1972,0.1874,0.1780,0.1690,0.1604,0.1524,0.1450,0.1382,0.1320,0.1264,0.1214,0.1168,0.1128,0.1091,0.1059,0.1030,0.1004,0.0981,0.0960,0.0941,0.0923,0.0907,0.0892,0.0877,0.0864,0.0850,0.0837,0.0824,0.0812,0.0799,0.0787,0.0775,0.0763,0.0752,0.0741,0.0731,0.0721,0.0712,0.0704,0.0696,0.0689,0.0682,0.0676,0.0671,0.0667,0.0663,0.0660,0.0657,0.0656,0.0656,
                                                                 0.1930,0.1834,0.1741,0.1652,0.1568,0.1489,0.1417,0.1351,0.1291,0.1236,0.1187,0.1144,0.1105,0.1070,0.1039,0.1012,0.0987,0.0965,0.0945,0.0926,0.0910,0.0894,0.0879,0.0865,0.0852,0.0839,0.0826,0.0814,0.0802,0.0789,0.0777,0.0766,0.0754,0.0743,0.0733,0.0723,0.0714,0.0705,0.0697,0.0689,0.0682,0.0676,0.0670,0.0664,0.0659,0.0655,0.0652,0.0649,0.0648,0.0648,
                                                                 0.1890,0.1795,0.1703,0.1616,0.1533,0.1457,0.1386,0.1322,0.1263,0.1211,0.1164,0.1121,0.1084,0.1051,0.1021,0.0995,0.0971,0.0950,0.0931,0.0913,0.0897,0.0882,0.0868,0.0855,0.0842,0.0829,0.0816,0.0804,0.0792,0.0780,0.0768,0.0757,0.0746,0.0735,0.0725,0.0715,0.0706,0.0698,0.0690,0.0682,0.0675,0.0669,0.0663,0.0658,0.0653,0.0648,0.0645,0.0642,0.0640,0.0639,
                                                                 0.1852,0.1758,0.1668,0.1582,0.1501,0.1426,0.1358,0.1295,0.1238,0.1187,0.1142,0.1101,0.1065,0.1033,0.1005,0.0980,0.0957,0.0937,0.0918,0.0902,0.0886,0.0872,0.0858,0.0845,0.0832,0.0819,0.0807,0.0795,0.0783,0.0771,0.0760,0.0749,0.0738,0.0728,0.0718,0.0708,0.0700,0.0691,0.0684,0.0676,0.0669,0.0663,0.0657,0.0651,0.0646,0.0641,0.0637,0.0634,0.0632,0.0632,
                                                                 0.1816,0.1723,0.1634,0.1550,0.1471,0.1398,0.1331,0.1270,0.1215,0.1166,0.1122,0.1083,0.1049,0.1018,0.0991,0.0966,0.0945,0.0925,0.0907,0.0891,0.0876,0.0862,0.0848,0.0835,0.0823,0.0811,0.0798,0.0787,0.0775,0.0763,0.0752,0.0741,0.0731,0.0720,0.0711,0.0702,0.0693,0.0685,0.0678,0.0670,0.0664,0.0657,0.0651,0.0645,0.0640,0.0635,0.0631,0.0627,0.0625,0.0624,
                                                                 0.1782,0.1691,0.1603,0.1521,0.1444,0.1372,0.1307,0.1248,0.1195,0.1147,0.1105,0.1067,0.1033,0.1004,0.0978,0.0954,0.0933,0.0914,0.0897,0.0881,0.0867,0.0853,0.0840,0.0827,0.0815,0.0802,0.0790,0.0779,0.0767,0.0756,0.0745,0.0734,0.0724,0.0714,0.0704,0.0696,0.0687,0.0679,0.0672,0.0665,0.0658,0.0651,0.0645,0.0639,0.0633,0.0628,0.0624,0.0620,0.0618,0.0617,
                                                                 0.1750,0.1660,0.1575,0.1494,0.1418,0.1349,0.1285,0.1228,0.1176,0.1130,0.1089,0.1052,0.1020,0.0991,0.0966,0.0943,0.0923,0.0905,0.0888,0.0873,0.0858,0.0845,0.0832,0.0819,0.0807,0.0795,0.0783,0.0771,0.0759,0.0748,0.0737,0.0727,0.0717,0.0707,0.0698,0.0690,0.0682,0.0674,0.0667,0.0659,0.0653,0.0646,0.0639,0.0633,0.0627,0.0622,0.0617,0.0614,0.0611,0.0610,
                                                                 0.1720,0.1632,0.1548,0.1469,0.1395,0.1327,0.1265,0.1209,0.1159,0.1114,0.1074,0.1039,0.1008,0.0980,0.0955,0.0934,0.0914,0.0896,0.0880,0.0865,0.0850,0.0837,0.0824,0.0811,0.0799,0.0787,0.0775,0.0764,0.0752,0.0741,0.0731,0.0721,0.0711,0.0702,0.0693,0.0684,0.0676,0.0669,0.0661,0.0654,0.0647,0.0640,0.0634,0.0627,0.0621,0.0616,0.0611,0.0607,0.0604,0.0603,
                                                                 0.1693,0.1606,0.1524,0.1446,0.1374,0.1308,0.1247,0.1193,0.1144,0.1100,0.1061,0.1027,0.0997,0.0970,0.0946,0.0925,0.0905,0.0888,0.0872,0.0857,0.0843,0.0830,0.0817,0.0804,0.0792,0.0780,0.0768,0.0757,0.0746,0.0735,0.0724,0.0714,0.0705,0.0696,0.0687,0.0679,0.0671,0.0664,0.0656,0.0649,0.0642,0.0635,0.0628,0.0622,0.0616,0.0610,0.0605,0.0601,0.0598,0.0597,
                                                                 0.1668,0.1583,0.1501,0.1425,0.1355,0.1290,0.1231,0.1178,0.1130,0.1087,0.1050,0.1016,0.0987,0.0961,0.0937,0.0916,0.0898,0.0881,0.0865,0.0850,0.0836,0.0823,0.0810,0.0797,0.0785,0.0773,0.0762,0.0750,0.0739,0.0729,0.0718,0.0709,0.0699,0.0691,0.0682,0.0674,0.0666,0.0659,0.0652,0.0644,0.0637,0.0630,0.0623,0.0616,0.0610,0.0604,0.0599,0.0595,0.0592,0.0590,
                                                                 0.1645,0.1561,0.1481,0.1406,0.1337,0.1273,0.1216,0.1164,0.1117,0.1076,0.1039,0.1007,0.0978,0.0952,0.0929,0.0909,0.0890,0.0874,0.0858,0.0843,0.0829,0.0816,0.0803,0.0791,0.0779,0.0767,0.0755,0.0744,0.0733,0.0723,0.0713,0.0703,0.0694,0.0685,0.0677,0.0669,0.0662,0.0654,0.0647,0.0639,0.0632,0.0625,0.0618,0.0611,0.0604,0.0598,0.0593,0.0589,0.0586,0.0585,
                                                                 0.1623,0.1540,0.1462,0.1389,0.1321,0.1258,0.1202,0.1151,0.1106,0.1065,0.1029,0.0998,0.0969,0.0944,0.0922,0.0902,0.0884,0.0867,0.0851,0.0837,0.0823,0.0810,0.0797,0.0784,0.0772,0.0760,0.0749,0.0738,0.0727,0.0717,0.0707,0.0698,0.0689,0.0680,0.0672,0.0665,0.0657,0.0650,0.0642,0.0635,0.0627,0.0620,0.0613,0.0606,0.0599,0.0593,0.0588,0.0583,0.0580,0.0579,
                                                                 0.1603,0.1522,0.1445,0.1372,0.1306,0.1245,0.1189,0.1140,0.1095,0.1056,0.1020,0.0989,0.0962,0.0937,0.0915,0.0895,0.0877,0.0861,0.0845,0.0831,0.0817,0.0803,0.0791,0.0778,0.0766,0.0754,0.0743,0.0732,0.0721,0.0711,0.0702,0.0693,0.0684,0.0676,0.0668,0.0660,0.0653,0.0645,0.0638,0.0630,0.0622,0.0615,0.0608,0.0600,0.0594,0.0587,0.0582,0.0578,0.0575,0.0574,
                                                                 0.1585,0.1505,0.1429,0.1358,0.1292,0.1232,0.1178,0.1129,0.1085,0.1047,0.1012,0.0982,0.0954,0.0930,0.0909,0.0889,0.0871,0.0855,0.0839,0.0825,0.0811,0.0797,0.0784,0.0772,0.0760,0.0748,0.0737,0.0726,0.0716,0.0706,0.0697,0.0688,0.0679,0.0671,0.0663,0.0656,0.0648,0.0641,0.0633,0.0625,0.0618,0.0610,0.0602,0.0595,0.0588,0.0582,0.0577,0.0573,0.0570,0.0569,
                                                                 0.1568,0.1489,0.1414,0.1344,0.1279,0.1220,0.1167,0.1119,0.1076,0.1038,0.1004,0.0974,0.0948,0.0924,0.0902,0.0883,0.0865,0.0849,0.0833,0.0819,0.0805,0.0791,0.0778,0.0766,0.0754,0.0742,0.0731,0.0720,0.0710,0.0701,0.0692,0.0683,0.0675,0.0667,0.0659,0.0651,0.0644,0.0636,0.0629,0.0621,0.0613,0.0605,0.0598,0.0590,0.0583,0.0577,0.0572,0.0568,0.0565,0.0564,
                                                                 0.1552,0.1474,0.1400,0.1331,0.1267,0.1209,0.1157,0.1110,0.1068,0.1030,0.0997,0.0968,0.0941,0.0918,0.0896,0.0877,0.0859,0.0843,0.0827,0.0813,0.0799,0.0785,0.0772,0.0760,0.0748,0.0737,0.0726,0.0715,0.0705,0.0696,0.0687,0.0679,0.0670,0.0663,0.0655,0.0647,0.0640,0.0632,0.0624,0.0616,0.0608,0.0600,0.0593,0.0585,0.0578,0.0572,0.0567,0.0563,0.0561,0.0559,
                                                                 0.1537,0.1460,0.1387,0.1319,0.1256,0.1199,0.1147,0.1101,0.1060,0.1023,0.0990,0.0961,0.0935,0.0912,0.0891,0.0871,0.0854,0.0837,0.0822,0.0807,0.0793,0.0780,0.0767,0.0754,0.0742,0.0731,0.0720,0.0710,0.0700,0.0691,0.0683,0.0674,0.0666,0.0658,0.0651,0.0643,0.0635,0.0628,0.0620,0.0612,0.0604,0.0596,0.0588,0.0580,0.0574,0.0568,0.0563,0.0559,0.0556,0.0555,
                                                                 0.1523,0.1447,0.1374,0.1307,0.1246,0.1189,0.1138,0.1093,0.1052,0.1016,0.0984,0.0955,0.0929,0.0906,0.0885,0.0866,0.0848,0.0832,0.0816,0.0801,0.0787,0.0774,0.0761,0.0748,0.0737,0.0726,0.0715,0.0705,0.0696,0.0687,0.0678,0.0670,0.0662,0.0654,0.0647,0.0639,0.0631,0.0623,0.0615,0.0607,0.0599,0.0591,0.0583,0.0576,0.0569,0.0563,0.0558,0.0555,0.0552,0.0551,
                                                                 0.1510,0.1434,0.1363,0.1296,0.1236,0.1180,0.1130,0.1085,0.1045,0.1009,0.0977,0.0949,0.0923,0.0900,0.0879,0.0860,0.0842,0.0826,0.0810,0.0795,0.0781,0.0768,0.0755,0.0743,0.0731,0.0720,0.0710,0.0700,0.0691,0.0682,0.0674,0.0666,0.0658,0.0650,0.0643,0.0635,0.0627,0.0619,0.0611,0.0602,0.0594,0.0586,0.0578,0.0571,0.0565,0.0559,0.0554,0.0551,0.0548,0.0547,
                                                                 0.1498,0.1422,0.1352,0.1286,0.1226,0.1171,0.1122,0.1077,0.1038,0.1002,0.0971,0.0943,0.0917,0.0895,0.0874,0.0855,0.0837,0.0820,0.0805,0.0790,0.0776,0.0762,0.0750,0.0737,0.0726,0.0715,0.0705,0.0696,0.0687,0.0678,0.0670,0.0662,0.0654,0.0646,0.0639,0.0631,0.0623,0.0615,0.0606,0.0598,0.0590,0.0582,0.0574,0.0567,0.0560,0.0555,0.0550,0.0547,0.0545,0.0544,
                                                                 0.1486,0.1411,0.1341,0.1276,0.1217,0.1163,0.1114,0.1070,0.1031,0.0996,0.0965,0.0937,0.0912,0.0889,0.0868,0.0849,0.0831,0.0815,0.0799,0.0784,0.0770,0.0757,0.0744,0.0732,0.0721,0.0710,0.0700,0.0691,0.0682,0.0674,0.0666,0.0658,0.0650,0.0642,0.0635,0.0627,0.0619,0.0610,0.0602,0.0593,0.0585,0.0577,0.0569,0.0562,0.0556,0.0551,0.0546,0.0543,0.0541,0.0540,
                                                                 0.1474,0.1400,0.1331,0.1267,0.1208,0.1155,0.1106,0.1063,0.1024,0.0990,0.0959,0.0931,0.0906,0.0883,0.0863,0.0843,0.0826,0.0809,0.0793,0.0778,0.0764,0.0751,0.0739,0.0727,0.0716,0.0705,0.0696,0.0687,0.0678,0.0670,0.0662,0.0654,0.0646,0.0638,0.0631,0.0622,0.0614,0.0606,0.0597,0.0589,0.0581,0.0573,0.0565,0.0558,0.0552,0.0547,0.0543,0.0540,0.0538,0.0537,
                                                                 0.1463,0.1390,0.1321,0.1257,0.1199,0.1147,0.1099,0.1056,0.1018,0.0984,0.0953,0.0925,0.0900,0.0878,0.0857,0.0838,0.0820,0.0803,0.0787,0.0773,0.0759,0.0746,0.0733,0.0722,0.0711,0.0701,0.0691,0.0682,0.0674,0.0666,0.0658,0.0650,0.0642,0.0635,0.0627,0.0618,0.0610,0.0601,0.0593,0.0584,0.0576,0.0568,0.0561,0.0554,0.0548,0.0543,0.0540,0.0537,0.0535,0.0534,
                                                                 0.1452,0.1379,0.1311,0.1248,0.1191,0.1139,0.1092,0.1049,0.1011,0.0977,0.0947,0.0919,0.0895,0.0872,0.0851,0.0832,0.0814,0.0797,0.0782,0.0767,0.0753,0.0740,0.0728,0.0717,0.0706,0.0696,0.0687,0.0678,0.0670,0.0662,0.0654,0.0646,0.0639,0.0631,0.0622,0.0614,0.0606,0.0597,0.0588,0.0580,0.0572,0.0564,0.0557,0.0550,0.0545,0.0540,0.0536,0.0534,0.0532,0.0532,
                                                                 0.1442,0.1369,0.1302,0.1240,0.1183,0.1131,0.1084,0.1043,0.1005,0.0971,0.0941,0.0914,0.0889,0.0866,0.0845,0.0826,0.0808,0.0792,0.0776,0.0761,0.0748,0.0735,0.0723,0.0712,0.0702,0.0692,0.0683,0.0674,0.0666,0.0658,0.0650,0.0643,0.0635,0.0627,0.0618,0.0610,0.0601,0.0593,0.0584,0.0576,0.0568,0.0560,0.0553,0.0547,0.0541,0.0537,0.0533,0.0531,0.0530,0.0529,
                                                                 0.1431,0.1360,0.1293,0.1231,0.1175,0.1124,0.1077,0.1036,0.0999,0.0965,0.0935,0.0908,0.0883,0.0860,0.0840,0.0820,0.0803,0.0786,0.0770,0.0756,0.0742,0.0730,0.0718,0.0707,0.0697,0.0688,0.0679,0.0670,0.0662,0.0654,0.0647,0.0639,0.0631,0.0623,0.0614,0.0606,0.0597,0.0588,0.0580,0.0571,0.0563,0.0556,0.0549,0.0543,0.0538,0.0534,0.0531,0.0528,0.0527,0.0527,
                                                                 0.1421,0.1350,0.1284,0.1223,0.1167,0.1116,0.1070,0.1029,0.0992,0.0959,0.0929,0.0902,0.0877,0.0854,0.0834,0.0815,0.0797,0.0780,0.0765,0.0750,0.0737,0.0725,0.0713,0.0703,0.0693,0.0684,0.0675,0.0667,0.0659,0.0651,0.0643,0.0635,0.0627,0.0619,0.0610,0.0601,0.0593,0.0584,0.0575,0.0567,0.0559,0.0552,0.0546,0.0540,0.0535,0.0531,0.0528,0.0526,0.0525,0.0525,
                                                                 0.1411,0.1341,0.1275,0.1214,0.1159,0.1109,0.1063,0.1022,0.0986,0.0953,0.0923,0.0896,0.0871,0.0849,0.0828,0.0809,0.0791,0.0774,0.0759,0.0745,0.0732,0.0720,0.0709,0.0698,0.0689,0.0680,0.0671,0.0663,0.0655,0.0647,0.0639,0.0631,0.0623,0.0614,0.0606,0.0597,0.0588,0.0580,0.0571,0.0563,0.0556,0.0549,0.0542,0.0537,0.0532,0.0529,0.0526,0.0524,0.0523,0.0522,
                                                                 0.1401,0.1331,0.1266,0.1206,0.1151,0.1101,0.1056,0.1016,0.0979,0.0946,0.0917,0.0890,0.0865,0.0843,0.0822,0.0803,0.0785,0.0769,0.0754,0.0740,0.0727,0.0715,0.0704,0.0694,0.0685,0.0676,0.0667,0.0659,0.0651,0.0643,0.0635,0.0627,0.0619,0.0610,0.0602,0.0593,0.0584,0.0575,0.0567,0.0559,0.0552,0.0545,0.0539,0.0534,0.0530,0.0526,0.0523,0.0522,0.0521,0.0521,
                                                                 0.1392,0.1322,0.1257,0.1198,0.1143,0.1094,0.1049,0.1009,0.0973,0.0940,0.0910,0.0883,0.0859,0.0836,0.0816,0.0797,0.0779,0.0763,0.0748,0.0735,0.0722,0.0710,0.0700,0.0690,0.0681,0.0672,0.0664,0.0655,0.0647,0.0640,0.0631,0.0623,0.0615,0.0606,0.0597,0.0589,0.0580,0.0571,0.0563,0.0555,0.0548,0.0542,0.0536,0.0531,0.0527,0.0524,0.0521,0.0520,0.0519,0.0519,
                                                                 0.1382,0.1313,0.1248,0.1189,0.1135,0.1086,0.1042,0.1002,0.0966,0.0933,0.0904,0.0877,0.0853,0.0830,0.0810,0.0791,0.0774,0.0758,0.0743,0.0730,0.0717,0.0706,0.0695,0.0686,0.0677,0.0668,0.0660,0.0652,0.0644,0.0636,0.0628,0.0619,0.0611,0.0602,0.0593,0.0584,0.0576,0.0567,0.0559,0.0552,0.0545,0.0539,0.0533,0.0529,0.0525,0.0522,0.0519,0.0518,0.0517,0.0517,
                                                                 0.1372,0.1304,0.1240,0.1181,0.1128,0.1079,0.1035,0.0995,0.0959,0.0927,0.0897,0.0871,0.0846,0.0824,0.0804,0.0785,0.0768,0.0752,0.0738,0.0725,0.0713,0.0701,0.0691,0.0682,0.0673,0.0664,0.0656,0.0648,0.0640,0.0632,0.0624,0.0615,0.0607,0.0598,0.0589,0.0580,0.0571,0.0563,0.0555,0.0548,0.0541,0.0536,0.0530,0.0526,0.0523,0.0520,0.0518,0.0516,0.0516,0.0516};
        
        return arrayCurvatureMax[y * arrayWidth + x];
    }
    case 10:
    {
        static const std::vector<double> arrayCurvatureMax =    {0.2448,0.2377,0.2304,0.2229,0.2156,0.2083,0.2014,0.1948,0.1887,0.1829,0.1776,0.1727,0.1682,0.1639,0.1600,0.1563,0.1528,0.1494,0.1462,0.1431,0.1402,0.1373,0.1344,0.1317,0.1290,0.1264,0.1238,0.1212,0.1188,0.1163,0.1140,0.1117,0.1094,0.1072,0.1050,0.1029,0.1009,0.0989,0.0970,0.0951,0.0933,0.0916,0.0899,0.0883,0.0867,0.0852,0.0838,0.0824,0.0811,0.0798,
                                                                 0.2424,0.2352,0.2277,0.2202,0.2127,0.2054,0.1984,0.1918,0.1856,0.1798,0.1744,0.1695,0.1649,0.1606,0.1566,0.1529,0.1493,0.1460,0.1427,0.1396,0.1367,0.1338,0.1309,0.1282,0.1255,0.1229,0.1204,0.1179,0.1154,0.1131,0.1107,0.1085,0.1063,0.1041,0.1020,0.1000,0.0981,0.0962,0.0943,0.0925,0.0908,0.0892,0.0876,0.0860,0.0846,0.0831,0.0818,0.0805,0.0792,0.0780,
                                                                 0.2399,0.2325,0.2250,0.2173,0.2098,0.2024,0.1953,0.1886,0.1824,0.1765,0.1711,0.1661,0.1615,0.1572,0.1532,0.1494,0.1459,0.1425,0.1393,0.1362,0.1332,0.1303,0.1275,0.1247,0.1221,0.1195,0.1170,0.1146,0.1122,0.1099,0.1076,0.1054,0.1033,0.1012,0.0992,0.0972,0.0954,0.0935,0.0918,0.0901,0.0885,0.0869,0.0854,0.0839,0.0825,0.0812,0.0799,0.0786,0.0775,0.0763,
                                                                 0.2373,0.2298,0.2221,0.2144,0.2067,0.1993,0.1922,0.1854,0.1791,0.1732,0.1678,0.1628,0.1581,0.1538,0.1498,0.1460,0.1424,0.1390,0.1358,0.1327,0.1297,0.1268,0.1240,0.1214,0.1187,0.1162,0.1137,0.1113,0.1090,0.1068,0.1046,0.1024,0.1004,0.0984,0.0964,0.0946,0.0928,0.0910,0.0894,0.0878,0.0862,0.0847,0.0833,0.0819,0.0806,0.0793,0.0781,0.0769,0.0758,0.0748,
                                                                 0.2345,0.2269,0.2191,0.2113,0.2036,0.1961,0.1889,0.1821,0.1758,0.1699,0.1644,0.1594,0.1547,0.1503,0.1463,0.1425,0.1389,0.1355,0.1323,0.1292,0.1263,0.1234,0.1207,0.1180,0.1155,0.1130,0.1106,0.1082,0.1060,0.1038,0.1016,0.0996,0.0976,0.0957,0.0938,0.0920,0.0903,0.0887,0.0871,0.0855,0.0841,0.0827,0.0813,0.0800,0.0788,0.0776,0.0765,0.0754,0.0743,0.0733,
                                                                 0.2316,0.2239,0.2161,0.2082,0.2004,0.1928,0.1856,0.1787,0.1724,0.1664,0.1610,0.1559,0.1512,0.1469,0.1428,0.1390,0.1355,0.1321,0.1289,0.1258,0.1229,0.1201,0.1174,0.1148,0.1123,0.1098,0.1075,0.1052,0.1030,0.1009,0.0988,0.0969,0.0949,0.0931,0.0913,0.0896,0.0880,0.0864,0.0849,0.0835,0.0821,0.0807,0.0795,0.0783,0.0771,0.0760,0.0749,0.0739,0.0729,0.0720,
                                                                 0.2287,0.2209,0.2129,0.2049,0.1971,0.1894,0.1822,0.1753,0.1689,0.1630,0.1575,0.1524,0.1478,0.1434,0.1394,0.1356,0.1321,0.1287,0.1255,0.1225,0.1196,0.1168,0.1142,0.1116,0.1092,0.1068,0.1045,0.1023,0.1002,0.0981,0.0962,0.0943,0.0924,0.0907,0.0890,0.0874,0.0858,0.0843,0.0829,0.0815,0.0802,0.0789,0.0777,0.0766,0.0755,0.0745,0.0735,0.0725,0.0716,0.0707,
                                                                 0.2256,0.2177,0.2096,0.2016,0.1937,0.1860,0.1787,0.1719,0.1654,0.1595,0.1540,0.1490,0.1443,0.1400,0.1360,0.1322,0.1287,0.1254,0.1222,0.1192,0.1164,0.1137,0.1111,0.1086,0.1062,0.1039,0.1016,0.0995,0.0975,0.0955,0.0936,0.0918,0.0900,0.0884,0.0867,0.0852,0.0837,0.0823,0.0810,0.0797,0.0784,0.0773,0.0761,0.0751,0.0740,0.0731,0.0721,0.0712,0.0704,0.0695,
                                                                 0.2224,0.2144,0.2063,0.1982,0.1902,0.1826,0.1752,0.1684,0.1620,0.1560,0.1506,0.1455,0.1409,0.1366,0.1326,0.1288,0.1254,0.1221,0.1190,0.1160,0.1132,0.1106,0.1080,0.1056,0.1033,0.1010,0.0989,0.0968,0.0949,0.0930,0.0912,0.0894,0.0878,0.0862,0.0846,0.0832,0.0818,0.0805,0.0792,0.0780,0.0768,0.0757,0.0746,0.0736,0.0727,0.0717,0.0709,0.0700,0.0692,0.0684,
                                                                 0.2192,0.2111,0.2029,0.1948,0.1868,0.1791,0.1717,0.1649,0.1585,0.1526,0.1471,0.1421,0.1375,0.1332,0.1292,0.1256,0.1221,0.1189,0.1158,0.1130,0.1102,0.1076,0.1051,0.1028,0.1005,0.0984,0.0963,0.0943,0.0924,0.0906,0.0889,0.0872,0.0856,0.0841,0.0827,0.0813,0.0800,0.0787,0.0775,0.0764,0.0753,0.0742,0.0732,0.0723,0.0714,0.0705,0.0697,0.0689,0.0681,0.0674,
                                                                 0.2159,0.2077,0.1995,0.1913,0.1833,0.1756,0.1682,0.1614,0.1550,0.1491,0.1437,0.1387,0.1341,0.1299,0.1260,0.1223,0.1190,0.1158,0.1128,0.1100,0.1073,0.1048,0.1024,0.1001,0.0979,0.0958,0.0938,0.0919,0.0901,0.0884,0.0867,0.0851,0.0836,0.0822,0.0808,0.0795,0.0783,0.0771,0.0759,0.0749,0.0738,0.0729,0.0719,0.0710,0.0702,0.0694,0.0686,0.0679,0.0671,0.0665,
                                                                 0.2125,0.2043,0.1960,0.1878,0.1798,0.1720,0.1647,0.1579,0.1516,0.1457,0.1403,0.1354,0.1308,0.1267,0.1228,0.1192,0.1159,0.1128,0.1099,0.1071,0.1045,0.1020,0.0997,0.0975,0.0954,0.0934,0.0915,0.0896,0.0879,0.0863,0.0847,0.0832,0.0817,0.0804,0.0791,0.0778,0.0767,0.0756,0.0745,0.0735,0.0725,0.0716,0.0707,0.0699,0.0691,0.0683,0.0676,0.0669,0.0662,0.0656,
                                                                 0.2091,0.2008,0.1925,0.1843,0.1762,0.1685,0.1613,0.1545,0.1482,0.1424,0.1370,0.1321,0.1277,0.1235,0.1197,0.1162,0.1129,0.1099,0.1070,0.1044,0.1018,0.0994,0.0972,0.0950,0.0930,0.0911,0.0892,0.0875,0.0858,0.0843,0.0828,0.0813,0.0800,0.0787,0.0775,0.0763,0.0752,0.0741,0.0731,0.0722,0.0713,0.0704,0.0696,0.0688,0.0681,0.0673,0.0666,0.0660,0.0654,0.0647,
                                                                 0.2057,0.1974,0.1890,0.1808,0.1728,0.1651,0.1578,0.1511,0.1448,0.1391,0.1338,0.1290,0.1246,0.1205,0.1168,0.1133,0.1101,0.1071,0.1043,0.1017,0.0993,0.0970,0.0948,0.0927,0.0908,0.0889,0.0872,0.0855,0.0839,0.0824,0.0810,0.0796,0.0783,0.0771,0.0760,0.0749,0.0738,0.0728,0.0719,0.0710,0.0701,0.0693,0.0685,0.0678,0.0671,0.0664,0.0658,0.0651,0.0645,0.0640,
                                                                 0.2022,0.1939,0.1855,0.1773,0.1693,0.1617,0.1545,0.1478,0.1416,0.1359,0.1307,0.1259,0.1216,0.1176,0.1139,0.1105,0.1074,0.1045,0.1018,0.0992,0.0969,0.0946,0.0925,0.0905,0.0887,0.0869,0.0852,0.0836,0.0821,0.0807,0.0793,0.0780,0.0768,0.0756,0.0746,0.0735,0.0725,0.0716,0.0707,0.0699,0.0690,0.0683,0.0676,0.0669,0.0662,0.0656,0.0649,0.0643,0.0638,0.0632,
                                                                 0.1988,0.1904,0.1821,0.1739,0.1659,0.1583,0.1512,0.1445,0.1384,0.1328,0.1276,0.1230,0.1187,0.1148,0.1112,0.1079,0.1048,0.1020,0.0993,0.0969,0.0946,0.0924,0.0904,0.0885,0.0867,0.0850,0.0834,0.0818,0.0804,0.0790,0.0778,0.0765,0.0754,0.0743,0.0732,0.0723,0.0713,0.0704,0.0696,0.0688,0.0681,0.0673,0.0666,0.0660,0.0654,0.0648,0.0642,0.0636,0.0631,0.0625,
                                                                 0.1954,0.1870,0.1787,0.1705,0.1626,0.1550,0.1479,0.1414,0.1353,0.1298,0.1247,0.1201,0.1159,0.1121,0.1086,0.1053,0.1024,0.0996,0.0970,0.0947,0.0924,0.0904,0.0884,0.0866,0.0848,0.0832,0.0817,0.0802,0.0788,0.0775,0.0763,0.0751,0.0741,0.0730,0.0720,0.0711,0.0702,0.0694,0.0686,0.0678,0.0671,0.0664,0.0658,0.0652,0.0646,0.0640,0.0634,0.0629,0.0624,0.0618,
                                                                 0.1920,0.1836,0.1753,0.1672,0.1593,0.1518,0.1448,0.1383,0.1323,0.1269,0.1219,0.1174,0.1133,0.1095,0.1061,0.1029,0.1000,0.0974,0.0949,0.0926,0.0904,0.0884,0.0865,0.0848,0.0831,0.0815,0.0801,0.0787,0.0774,0.0761,0.0750,0.0739,0.0728,0.0718,0.0709,0.0700,0.0692,0.0684,0.0677,0.0669,0.0663,0.0656,0.0650,0.0644,0.0638,0.0633,0.0627,0.0622,0.0617,0.0612,
                                                                 0.1886,0.1803,0.1720,0.1639,0.1561,0.1487,0.1418,0.1354,0.1295,0.1241,0.1192,0.1148,0.1108,0.1071,0.1037,0.1007,0.0979,0.0953,0.0929,0.0906,0.0885,0.0866,0.0848,0.0831,0.0815,0.0800,0.0786,0.0773,0.0760,0.0748,0.0737,0.0727,0.0717,0.0707,0.0699,0.0690,0.0682,0.0675,0.0668,0.0661,0.0654,0.0648,0.0642,0.0637,0.0631,0.0626,0.0621,0.0616,0.0611,0.0606,
                                                                 0.1853,0.1771,0.1688,0.1608,0.1530,0.1457,0.1389,0.1325,0.1268,0.1215,0.1167,0.1123,0.1084,0.1048,0.1015,0.0986,0.0958,0.0933,0.0910,0.0888,0.0868,0.0849,0.0832,0.0815,0.0800,0.0786,0.0772,0.0759,0.0747,0.0736,0.0726,0.0716,0.0706,0.0697,0.0689,0.0681,0.0673,0.0666,0.0660,0.0653,0.0647,0.0641,0.0635,0.0630,0.0624,0.0619,0.0614,0.0609,0.0604,0.0599,
                                                                 0.1821,0.1739,0.1657,0.1577,0.1501,0.1428,0.1361,0.1298,0.1242,0.1190,0.1143,0.1100,0.1062,0.1027,0.0995,0.0966,0.0939,0.0914,0.0892,0.0871,0.0851,0.0833,0.0817,0.0801,0.0786,0.0772,0.0759,0.0747,0.0736,0.0725,0.0715,0.0705,0.0696,0.0688,0.0680,0.0672,0.0665,0.0658,0.0652,0.0646,0.0640,0.0634,0.0628,0.0623,0.0618,0.0613,0.0608,0.0603,0.0598,0.0593,
                                                                 0.1790,0.1708,0.1627,0.1548,0.1472,0.1401,0.1334,0.1273,0.1217,0.1166,0.1120,0.1078,0.1041,0.1006,0.0975,0.0947,0.0921,0.0897,0.0875,0.0855,0.0836,0.0819,0.0803,0.0787,0.0773,0.0760,0.0748,0.0736,0.0725,0.0715,0.0705,0.0696,0.0687,0.0679,0.0672,0.0664,0.0657,0.0651,0.0645,0.0639,0.0633,0.0627,0.0622,0.0617,0.0612,0.0607,0.0602,0.0597,0.0592,0.0587,
                                                                 0.1760,0.1679,0.1598,0.1520,0.1445,0.1374,0.1309,0.1249,0.1194,0.1144,0.1099,0.1058,0.1021,0.0988,0.0957,0.0930,0.0905,0.0881,0.0860,0.0841,0.0822,0.0805,0.0790,0.0775,0.0761,0.0749,0.0737,0.0726,0.0715,0.0705,0.0696,0.0687,0.0679,0.0671,0.0664,0.0657,0.0650,0.0644,0.0638,0.0632,0.0626,0.0621,0.0616,0.0611,0.0606,0.0601,0.0596,0.0591,0.0586,0.0581,
                                                                 0.1731,0.1650,0.1570,0.1493,0.1419,0.1349,0.1285,0.1226,0.1172,0.1123,0.1079,0.1039,0.1003,0.0970,0.0941,0.0914,0.0889,0.0867,0.0846,0.0827,0.0809,0.0793,0.0778,0.0764,0.0751,0.0738,0.0727,0.0716,0.0706,0.0696,0.0687,0.0679,0.0671,0.0663,0.0656,0.0650,0.0643,0.0637,0.0631,0.0625,0.0620,0.0615,0.0609,0.0604,0.0599,0.0595,0.0590,0.0585,0.0580,0.0575,
                                                                 0.1703,0.1623,0.1544,0.1467,0.1394,0.1326,0.1262,0.1204,0.1151,0.1103,0.1060,0.1021,0.0986,0.0954,0.0925,0.0899,0.0875,0.0853,0.0833,0.0815,0.0798,0.0782,0.0767,0.0753,0.0741,0.0729,0.0718,0.0707,0.0697,0.0688,0.0679,0.0671,0.0664,0.0656,0.0649,0.0643,0.0636,0.0630,0.0625,0.0619,0.0614,0.0609,0.0603,0.0598,0.0593,0.0589,0.0584,0.0579,0.0574,0.0569,
                                                                 0.1676,0.1596,0.1518,0.1443,0.1371,0.1303,0.1241,0.1184,0.1132,0.1085,0.1043,0.1005,0.0970,0.0939,0.0911,0.0886,0.0862,0.0841,0.0821,0.0803,0.0787,0.0771,0.0757,0.0744,0.0731,0.0720,0.0709,0.0699,0.0689,0.0680,0.0672,0.0664,0.0657,0.0649,0.0643,0.0636,0.0630,0.0624,0.0619,0.0613,0.0608,0.0602,0.0597,0.0592,0.0587,0.0582,0.0577,0.0573,0.0567,0.0562,
                                                                 0.1650,0.1571,0.1494,0.1420,0.1349,0.1282,0.1221,0.1165,0.1114,0.1068,0.1027,0.0989,0.0956,0.0925,0.0898,0.0873,0.0850,0.0830,0.0811,0.0793,0.0777,0.0762,0.0748,0.0735,0.0723,0.0712,0.0701,0.0691,0.0682,0.0673,0.0665,0.0657,0.0650,0.0643,0.0636,0.0630,0.0624,0.0618,0.0612,0.0607,0.0602,0.0596,0.0591,0.0586,0.0581,0.0576,0.0571,0.0566,0.0561,0.0556,
                                                                 0.1625,0.1548,0.1471,0.1398,0.1328,0.1263,0.1202,0.1147,0.1097,0.1052,0.1012,0.0975,0.0942,0.0913,0.0886,0.0862,0.0840,0.0819,0.0801,0.0784,0.0768,0.0753,0.0740,0.0727,0.0715,0.0704,0.0694,0.0684,0.0675,0.0667,0.0659,0.0651,0.0644,0.0637,0.0630,0.0624,0.0618,0.0612,0.0606,0.0601,0.0596,0.0590,0.0585,0.0580,0.0575,0.0570,0.0565,0.0560,0.0555,0.0550,
                                                                 0.1602,0.1525,0.1450,0.1377,0.1309,0.1244,0.1185,0.1131,0.1082,0.1038,0.0998,0.0963,0.0930,0.0901,0.0875,0.0851,0.0830,0.0810,0.0792,0.0775,0.0759,0.0745,0.0732,0.0719,0.0708,0.0697,0.0687,0.0678,0.0669,0.0660,0.0652,0.0645,0.0638,0.0631,0.0624,0.0618,0.0612,0.0606,0.0601,0.0595,0.0590,0.0585,0.0579,0.0574,0.0569,0.0564,0.0559,0.0554,0.0549,0.0544,
                                                                 0.1580,0.1504,0.1430,0.1358,0.1290,0.1227,0.1169,0.1116,0.1068,0.1025,0.0986,0.0951,0.0919,0.0891,0.0865,0.0842,0.0821,0.0801,0.0783,0.0767,0.0752,0.0738,0.0725,0.0713,0.0701,0.0691,0.0681,0.0671,0.0662,0.0654,0.0646,0.0639,0.0632,0.0625,0.0618,0.0612,0.0606,0.0600,0.0595,0.0589,0.0584,0.0579,0.0573,0.0568,0.0563,0.0558,0.0553,0.0548,0.0543,0.0537,
                                                                 0.1560,0.1485,0.1411,0.1340,0.1274,0.1211,0.1154,0.1102,0.1055,0.1013,0.0974,0.0940,0.0909,0.0881,0.0856,0.0833,0.0813,0.0793,0.0776,0.0760,0.0745,0.0731,0.0718,0.0706,0.0695,0.0684,0.0675,0.0665,0.0657,0.0648,0.0640,0.0633,0.0626,0.0619,0.0613,0.0606,0.0600,0.0594,0.0589,0.0583,0.0578,0.0573,0.0567,0.0562,0.0557,0.0552,0.0547,0.0542,0.0537,0.0532,
                                                                 0.1540,0.1466,0.1393,0.1324,0.1258,0.1197,0.1141,0.1089,0.1043,0.1002,0.0964,0.0930,0.0900,0.0873,0.0848,0.0826,0.0805,0.0786,0.0769,0.0753,0.0738,0.0725,0.0712,0.0700,0.0689,0.0679,0.0669,0.0660,0.0651,0.0643,0.0635,0.0627,0.0620,0.0613,0.0607,0.0601,0.0594,0.0589,0.0583,0.0577,0.0572,0.0567,0.0561,0.0556,0.0551,0.0546,0.0541,0.0536,0.0531,0.0526,
                                                                 0.1522,0.1449,0.1377,0.1308,0.1243,0.1183,0.1128,0.1078,0.1032,0.0991,0.0955,0.0922,0.0892,0.0865,0.0841,0.0819,0.0798,0.0780,0.0763,0.0747,0.0733,0.0719,0.0706,0.0695,0.0684,0.0673,0.0663,0.0654,0.0645,0.0637,0.0629,0.0622,0.0615,0.0608,0.0601,0.0595,0.0589,0.0583,0.0577,0.0572,0.0566,0.0561,0.0556,0.0551,0.0546,0.0541,0.0536,0.0531,0.0526,0.0521,
                                                                 0.1505,0.1432,0.1362,0.1294,0.1230,0.1171,0.1116,0.1067,0.1022,0.0982,0.0946,0.0914,0.0885,0.0858,0.0834,0.0812,0.0792,0.0774,0.0757,0.0742,0.0727,0.0714,0.0701,0.0689,0.0678,0.0668,0.0658,0.0649,0.0640,0.0632,0.0624,0.0616,0.0609,0.0602,0.0596,0.0589,0.0583,0.0577,0.0571,0.0566,0.0560,0.0555,0.0550,0.0545,0.0540,0.0535,0.0530,0.0526,0.0521,0.0516,
                                                                 0.1489,0.1417,0.1347,0.1280,0.1218,0.1159,0.1106,0.1057,0.1014,0.0974,0.0939,0.0907,0.0878,0.0852,0.0828,0.0807,0.0787,0.0769,0.0752,0.0736,0.0722,0.0709,0.0696,0.0684,0.0673,0.0663,0.0653,0.0644,0.0635,0.0626,0.0618,0.0611,0.0604,0.0597,0.0590,0.0584,0.0577,0.0572,0.0566,0.0560,0.0555,0.0550,0.0545,0.0540,0.0535,0.0530,0.0525,0.0521,0.0516,0.0511,
                                                                 0.1474,0.1403,0.1334,0.1268,0.1206,0.1149,0.1096,0.1048,0.1005,0.0967,0.0932,0.0900,0.0872,0.0846,0.0823,0.0801,0.0782,0.0764,0.0747,0.0732,0.0717,0.0704,0.0691,0.0679,0.0668,0.0658,0.0648,0.0639,0.0630,0.0621,0.0613,0.0606,0.0598,0.0591,0.0585,0.0578,0.0572,0.0566,0.0561,0.0555,0.0550,0.0545,0.0540,0.0535,0.0530,0.0526,0.0521,0.0516,0.0512,0.0507,
                                                                 0.1460,0.1390,0.1322,0.1257,0.1196,0.1139,0.1087,0.1040,0.0998,0.0960,0.0925,0.0894,0.0866,0.0841,0.0818,0.0797,0.0777,0.0759,0.0743,0.0727,0.0713,0.0699,0.0687,0.0675,0.0664,0.0653,0.0643,0.0634,0.0625,0.0616,0.0608,0.0600,0.0593,0.0586,0.0579,0.0573,0.0567,0.0561,0.0555,0.0550,0.0545,0.0540,0.0535,0.0531,0.0526,0.0521,0.0517,0.0512,0.0508,0.0503,
                                                                 0.1447,0.1377,0.1310,0.1246,0.1186,0.1130,0.1079,0.1033,0.0991,0.0954,0.0920,0.0889,0.0861,0.0836,0.0813,0.0792,0.0773,0.0755,0.0738,0.0723,0.0708,0.0695,0.0682,0.0670,0.0659,0.0648,0.0638,0.0629,0.0620,0.0611,0.0603,0.0595,0.0588,0.0581,0.0574,0.0568,0.0562,0.0556,0.0551,0.0546,0.0541,0.0536,0.0531,0.0526,0.0522,0.0518,0.0513,0.0509,0.0504,0.0500,
                                                                 0.1434,0.1366,0.1300,0.1236,0.1177,0.1122,0.1072,0.1026,0.0985,0.0948,0.0914,0.0884,0.0857,0.0832,0.0809,0.0788,0.0769,0.0751,0.0734,0.0719,0.0704,0.0691,0.0678,0.0666,0.0654,0.0644,0.0634,0.0624,0.0615,0.0606,0.0598,0.0590,0.0583,0.0576,0.0570,0.0563,0.0557,0.0552,0.0547,0.0541,0.0536,0.0532,0.0527,0.0523,0.0518,0.0514,0.0510,0.0506,0.0501,0.0497,
                                                                 0.1423,0.1355,0.1290,0.1227,0.1169,0.1114,0.1065,0.1020,0.0979,0.0943,0.0910,0.0880,0.0853,0.0828,0.0805,0.0784,0.0765,0.0747,0.0730,0.0715,0.0700,0.0687,0.0674,0.0662,0.0650,0.0639,0.0629,0.0619,0.0610,0.0602,0.0594,0.0586,0.0579,0.0572,0.0565,0.0559,0.0553,0.0548,0.0543,0.0538,0.0533,0.0528,0.0524,0.0520,0.0515,0.0511,0.0507,0.0503,0.0499,0.0494,
                                                                 0.1412,0.1345,0.1280,0.1219,0.1161,0.1107,0.1059,0.1014,0.0974,0.0938,0.0905,0.0876,0.0849,0.0824,0.0801,0.0781,0.0761,0.0743,0.0727,0.0711,0.0696,0.0683,0.0670,0.0657,0.0646,0.0635,0.0625,0.0615,0.0606,0.0597,0.0589,0.0582,0.0575,0.0568,0.0561,0.0555,0.0550,0.0544,0.0539,0.0534,0.0530,0.0525,0.0521,0.0517,0.0512,0.0508,0.0504,0.0500,0.0496,0.0492,
                                                                 0.1402,0.1336,0.1272,0.1210,0.1153,0.1101,0.1053,0.1009,0.0969,0.0933,0.0901,0.0872,0.0845,0.0820,0.0798,0.0777,0.0758,0.0740,0.0723,0.0707,0.0693,0.0679,0.0666,0.0654,0.0642,0.0631,0.0621,0.0611,0.0602,0.0594,0.0585,0.0578,0.0571,0.0564,0.0558,0.0552,0.0546,0.0541,0.0536,0.0531,0.0527,0.0523,0.0518,0.0514,0.0510,0.0506,0.0502,0.0498,0.0494,0.0490,
                                                                 0.1393,0.1327,0.1263,0.1203,0.1147,0.1095,0.1047,0.1004,0.0965,0.0929,0.0897,0.0868,0.0841,0.0817,0.0795,0.0774,0.0755,0.0737,0.0720,0.0704,0.0689,0.0675,0.0662,0.0650,0.0638,0.0627,0.0617,0.0608,0.0599,0.0590,0.0582,0.0575,0.0568,0.0561,0.0555,0.0549,0.0544,0.0538,0.0534,0.0529,0.0525,0.0520,0.0516,0.0512,0.0508,0.0504,0.0500,0.0496,0.0492,0.0488,
                                                                 0.1383,0.1318,0.1255,0.1196,0.1140,0.1089,0.1042,0.0999,0.0960,0.0925,0.0893,0.0864,0.0838,0.0814,0.0791,0.0771,0.0751,0.0733,0.0717,0.0701,0.0686,0.0672,0.0659,0.0647,0.0635,0.0624,0.0614,0.0604,0.0595,0.0587,0.0579,0.0572,0.0565,0.0558,0.0552,0.0547,0.0541,0.0536,0.0531,0.0527,0.0523,0.0518,0.0514,0.0510,0.0506,0.0502,0.0498,0.0494,0.0490,0.0486,
                                                                 0.1375,0.1310,0.1248,0.1189,0.1134,0.1083,0.1037,0.0994,0.0956,0.0921,0.0890,0.0861,0.0835,0.0811,0.0788,0.0768,0.0748,0.0731,0.0714,0.0698,0.0683,0.0669,0.0656,0.0644,0.0632,0.0621,0.0611,0.0602,0.0593,0.0584,0.0577,0.0569,0.0563,0.0556,0.0550,0.0545,0.0539,0.0534,0.0530,0.0525,0.0521,0.0517,0.0513,0.0508,0.0504,0.0500,0.0496,0.0492,0.0488,0.0484,
                                                                 0.1366,0.1302,0.1241,0.1182,0.1128,0.1078,0.1032,0.0990,0.0952,0.0918,0.0886,0.0858,0.0832,0.0808,0.0786,0.0765,0.0746,0.0728,0.0711,0.0695,0.0681,0.0667,0.0654,0.0641,0.0630,0.0619,0.0609,0.0600,0.0591,0.0582,0.0575,0.0568,0.0561,0.0554,0.0549,0.0543,0.0538,0.0533,0.0528,0.0524,0.0519,0.0515,0.0511,0.0507,0.0503,0.0499,0.0495,0.0491,0.0487,0.0483,
                                                                 0.1358,0.1295,0.1234,0.1176,0.1122,0.1072,0.1027,0.0986,0.0948,0.0914,0.0883,0.0855,0.0829,0.0805,0.0783,0.0763,0.0744,0.0726,0.0709,0.0693,0.0678,0.0665,0.0652,0.0639,0.0628,0.0617,0.0607,0.0598,0.0589,0.0581,0.0573,0.0566,0.0559,0.0553,0.0547,0.0542,0.0537,0.0532,0.0527,0.0523,0.0518,0.0514,0.0510,0.0506,0.0502,0.0498,0.0494,0.0490,0.0486,0.0482,
                                                                 0.1350,0.1287,0.1227,0.1169,0.1116,0.1067,0.1022,0.0981,0.0944,0.0911,0.0880,0.0852,0.0827,0.0803,0.0781,0.0761,0.0742,0.0724,0.0707,0.0691,0.0677,0.0663,0.0650,0.0638,0.0627,0.0616,0.0606,0.0597,0.0588,0.0580,0.0572,0.0565,0.0558,0.0552,0.0546,0.0541,0.0536,0.0531,0.0526,0.0522,0.0518,0.0513,0.0509,0.0505,0.0501,0.0497,0.0493,0.0489,0.0485,0.0481,
                                                                 0.1342,0.1280,0.1220,0.1163,0.1111,0.1062,0.1018,0.0978,0.0941,0.0908,0.0877,0.0850,0.0824,0.0801,0.0779,0.0759,0.0740,0.0722,0.0706,0.0690,0.0676,0.0662,0.0649,0.0637,0.0626,0.0615,0.0605,0.0596,0.0587,0.0579,0.0572,0.0565,0.0558,0.0552,0.0546,0.0540,0.0535,0.0530,0.0526,0.0521,0.0517,0.0513,0.0508,0.0504,0.0500,0.0496,0.0492,0.0488,0.0484,0.0480,
                                                                 0.1335,0.1273,0.1213,0.1157,0.1105,0.1057,0.1013,0.0974,0.0938,0.0905,0.0875,0.0848,0.0822,0.0799,0.0778,0.0758,0.0739,0.0721,0.0705,0.0690,0.0675,0.0661,0.0649,0.0637,0.0625,0.0615,0.0605,0.0596,0.0587,0.0579,0.0571,0.0564,0.0558,0.0552,0.0546,0.0540,0.0535,0.0530,0.0525,0.0521,0.0516,0.0512,0.0508,0.0504,0.0499,0.0495,0.0491,0.0487,0.0483,0.0479};
        
        return arrayCurvatureMax[y * arrayWidth + x];
    }
    case 11:
    {
        static const std::vector<double> arrayCurvatureMax =    {0.2314,0.2225,0.2145,0.2073,0.2007,0.1948,0.1894,0.1845,0.1801,0.1759,0.1721,0.1685,0.1651,0.1617,0.1585,0.1553,0.1522,0.1490,0.1459,0.1427,0.1395,0.1364,0.1333,0.1302,0.1271,0.1242,0.1213,0.1185,0.1159,0.1134,0.1110,0.1087,0.1066,0.1046,0.1027,0.1010,0.0993,0.0977,0.0962,0.0947,0.0933,0.0919,0.0905,0.0892,0.0878,0.0866,0.0853,0.0840,0.0828,0.0816,
                                                                 0.2286,0.2197,0.2116,0.2043,0.1977,0.1917,0.1863,0.1813,0.1768,0.1726,0.1686,0.1649,0.1613,0.1579,0.1545,0.1512,0.1480,0.1447,0.1415,0.1382,0.1350,0.1318,0.1287,0.1256,0.1226,0.1197,0.1169,0.1142,0.1116,0.1092,0.1070,0.1048,0.1028,0.1009,0.0992,0.0975,0.0959,0.0944,0.0929,0.0915,0.0901,0.0888,0.0874,0.0861,0.0848,0.0836,0.0823,0.0811,0.0799,0.0788,
                                                                 0.2257,0.2168,0.2087,0.2014,0.1947,0.1887,0.1831,0.1781,0.1734,0.1691,0.1651,0.1612,0.1576,0.1540,0.1505,0.1471,0.1437,0.1404,0.1371,0.1338,0.1305,0.1273,0.1242,0.1211,0.1182,0.1153,0.1126,0.1100,0.1076,0.1053,0.1031,0.1011,0.0992,0.0974,0.0957,0.0942,0.0926,0.0912,0.0898,0.0884,0.0871,0.0858,0.0845,0.0832,0.0820,0.0807,0.0795,0.0783,0.0772,0.0760,
                                                                 0.2229,0.2140,0.2058,0.1984,0.1917,0.1856,0.1800,0.1748,0.1701,0.1657,0.1615,0.1576,0.1538,0.1501,0.1465,0.1430,0.1395,0.1361,0.1327,0.1294,0.1261,0.1229,0.1198,0.1168,0.1139,0.1111,0.1085,0.1061,0.1037,0.1016,0.0995,0.0976,0.0958,0.0941,0.0925,0.0910,0.0896,0.0882,0.0868,0.0855,0.0842,0.0830,0.0817,0.0805,0.0792,0.0780,0.0768,0.0757,0.0746,0.0735,
                                                                 0.2202,0.2112,0.2030,0.1955,0.1887,0.1825,0.1768,0.1716,0.1668,0.1622,0.1580,0.1539,0.1500,0.1462,0.1425,0.1389,0.1353,0.1318,0.1284,0.1251,0.1218,0.1186,0.1156,0.1126,0.1098,0.1072,0.1047,0.1023,0.1001,0.0980,0.0961,0.0943,0.0926,0.0910,0.0895,0.0881,0.0867,0.0854,0.0841,0.0828,0.0816,0.0803,0.0791,0.0779,0.0767,0.0755,0.0743,0.0732,0.0721,0.0711,
                                                                 0.2175,0.2085,0.2002,0.1927,0.1858,0.1795,0.1738,0.1684,0.1635,0.1588,0.1545,0.1503,0.1463,0.1424,0.1386,0.1349,0.1313,0.1277,0.1243,0.1209,0.1177,0.1146,0.1116,0.1087,0.1060,0.1035,0.1011,0.0988,0.0968,0.0948,0.0930,0.0913,0.0897,0.0882,0.0868,0.0854,0.0841,0.0828,0.0815,0.0803,0.0791,0.0779,0.0767,0.0755,0.0743,0.0732,0.0720,0.0709,0.0699,0.0688,
                                                                 0.2150,0.2059,0.1976,0.1900,0.1830,0.1766,0.1708,0.1653,0.1603,0.1555,0.1510,0.1467,0.1426,0.1386,0.1347,0.1310,0.1273,0.1238,0.1203,0.1170,0.1138,0.1107,0.1078,0.1050,0.1024,0.1000,0.0977,0.0956,0.0937,0.0918,0.0901,0.0885,0.0870,0.0856,0.0842,0.0829,0.0817,0.0804,0.0792,0.0780,0.0768,0.0757,0.0745,0.0733,0.0722,0.0710,0.0699,0.0688,0.0678,0.0668,
                                                                 0.2125,0.2034,0.1950,0.1873,0.1803,0.1738,0.1678,0.1623,0.1571,0.1523,0.1477,0.1433,0.1390,0.1350,0.1310,0.1272,0.1235,0.1200,0.1165,0.1132,0.1101,0.1071,0.1043,0.1016,0.0992,0.0969,0.0947,0.0927,0.0909,0.0891,0.0875,0.0860,0.0846,0.0832,0.0819,0.0807,0.0795,0.0783,0.0771,0.0759,0.0748,0.0736,0.0725,0.0713,0.0702,0.0691,0.0680,0.0669,0.0659,0.0650,
                                                                 0.2102,0.2010,0.1925,0.1848,0.1776,0.1711,0.1650,0.1594,0.1541,0.1491,0.1444,0.1399,0.1356,0.1314,0.1274,0.1236,0.1199,0.1163,0.1129,0.1097,0.1066,0.1038,0.1010,0.0985,0.0962,0.0940,0.0919,0.0901,0.0883,0.0867,0.0852,0.0837,0.0824,0.0811,0.0799,0.0787,0.0775,0.0763,0.0752,0.0740,0.0729,0.0718,0.0706,0.0695,0.0684,0.0673,0.0662,0.0652,0.0642,0.0633,
                                                                 0.2080,0.1987,0.1902,0.1823,0.1751,0.1684,0.1622,0.1565,0.1511,0.1460,0.1412,0.1366,0.1322,0.1280,0.1240,0.1201,0.1164,0.1129,0.1096,0.1064,0.1034,0.1007,0.0981,0.0956,0.0934,0.0913,0.0894,0.0877,0.0860,0.0845,0.0831,0.0817,0.0804,0.0792,0.0780,0.0768,0.0757,0.0746,0.0734,0.0723,0.0712,0.0701,0.0690,0.0679,0.0668,0.0657,0.0647,0.0637,0.0627,0.0618,
                                                                 0.2059,0.1965,0.1879,0.1800,0.1726,0.1658,0.1595,0.1537,0.1482,0.1430,0.1381,0.1334,0.1290,0.1247,0.1207,0.1168,0.1132,0.1097,0.1064,0.1034,0.1005,0.0978,0.0953,0.0930,0.0909,0.0890,0.0872,0.0855,0.0840,0.0825,0.0812,0.0799,0.0787,0.0775,0.0763,0.0752,0.0741,0.0730,0.0719,0.0708,0.0697,0.0685,0.0674,0.0664,0.0653,0.0642,0.0632,0.0622,0.0613,0.0604,
                                                                 0.2039,0.1944,0.1857,0.1776,0.1702,0.1633,0.1569,0.1509,0.1453,0.1400,0.1351,0.1303,0.1258,0.1216,0.1175,0.1137,0.1101,0.1067,0.1035,0.1005,0.0978,0.0952,0.0929,0.0907,0.0887,0.0868,0.0851,0.0836,0.0821,0.0808,0.0795,0.0782,0.0771,0.0759,0.0748,0.0737,0.0726,0.0715,0.0704,0.0693,0.0682,0.0672,0.0661,0.0650,0.0639,0.0629,0.0619,0.0610,0.0600,0.0592,
                                                                 0.2019,0.1924,0.1835,0.1754,0.1678,0.1608,0.1543,0.1482,0.1425,0.1372,0.1321,0.1273,0.1228,0.1186,0.1145,0.1107,0.1072,0.1039,0.1008,0.0979,0.0953,0.0928,0.0906,0.0885,0.0867,0.0849,0.0833,0.0818,0.0805,0.0792,0.0779,0.0768,0.0756,0.0745,0.0734,0.0724,0.0713,0.0702,0.0691,0.0681,0.0670,0.0659,0.0648,0.0638,0.0627,0.0617,0.0607,0.0598,0.0589,0.0581,
                                                                 0.2000,0.1903,0.1814,0.1731,0.1655,0.1583,0.1517,0.1455,0.1398,0.1343,0.1292,0.1244,0.1199,0.1157,0.1117,0.1080,0.1045,0.1013,0.0983,0.0955,0.0930,0.0907,0.0886,0.0866,0.0848,0.0832,0.0817,0.0803,0.0790,0.0777,0.0766,0.0754,0.0743,0.0733,0.0722,0.0711,0.0701,0.0690,0.0679,0.0669,0.0658,0.0647,0.0637,0.0626,0.0616,0.0606,0.0596,0.0587,0.0578,0.0570,
                                                                 0.1981,0.1883,0.1793,0.1709,0.1631,0.1559,0.1492,0.1429,0.1370,0.1316,0.1264,0.1216,0.1171,0.1129,0.1090,0.1053,0.1019,0.0988,0.0959,0.0933,0.0909,0.0887,0.0867,0.0849,0.0832,0.0816,0.0802,0.0789,0.0776,0.0764,0.0753,0.0742,0.0732,0.0721,0.0710,0.0700,0.0689,0.0679,0.0668,0.0658,0.0647,0.0636,0.0626,0.0615,0.0605,0.0596,0.0586,0.0577,0.0569,0.0561,
                                                                 0.1962,0.1863,0.1771,0.1686,0.1608,0.1534,0.1466,0.1403,0.1344,0.1289,0.1237,0.1189,0.1144,0.1103,0.1064,0.1028,0.0996,0.0966,0.0938,0.0913,0.0890,0.0869,0.0850,0.0833,0.0817,0.0802,0.0789,0.0776,0.0764,0.0753,0.0742,0.0731,0.0721,0.0710,0.0700,0.0689,0.0679,0.0668,0.0658,0.0647,0.0637,0.0626,0.0616,0.0605,0.0595,0.0586,0.0577,0.0568,0.0560,0.0552,
                                                                 0.1942,0.1843,0.1750,0.1664,0.1584,0.1510,0.1441,0.1377,0.1317,0.1262,0.1211,0.1163,0.1118,0.1078,0.1040,0.1005,0.0974,0.0945,0.0918,0.0894,0.0872,0.0853,0.0835,0.0818,0.0803,0.0789,0.0776,0.0764,0.0753,0.0742,0.0731,0.0721,0.0711,0.0700,0.0690,0.0680,0.0669,0.0659,0.0648,0.0637,0.0627,0.0616,0.0606,0.0596,0.0586,0.0577,0.0568,0.0559,0.0551,0.0544,
                                                                 0.1923,0.1822,0.1728,0.1641,0.1560,0.1485,0.1416,0.1351,0.1291,0.1236,0.1185,0.1137,0.1094,0.1054,0.1017,0.0983,0.0953,0.0925,0.0900,0.0877,0.0856,0.0838,0.0821,0.0805,0.0791,0.0778,0.0765,0.0754,0.0742,0.0732,0.0721,0.0711,0.0701,0.0691,0.0681,0.0670,0.0660,0.0649,0.0639,0.0628,0.0618,0.0607,0.0597,0.0587,0.0577,0.0568,0.0559,0.0551,0.0543,0.0536,
                                                                 0.1903,0.1801,0.1706,0.1618,0.1536,0.1461,0.1391,0.1326,0.1266,0.1211,0.1160,0.1113,0.1070,0.1031,0.0996,0.0963,0.0934,0.0907,0.0883,0.0861,0.0842,0.0824,0.0808,0.0793,0.0779,0.0767,0.0755,0.0744,0.0733,0.0722,0.0712,0.0702,0.0692,0.0682,0.0672,0.0661,0.0651,0.0640,0.0630,0.0619,0.0609,0.0598,0.0588,0.0579,0.0569,0.0560,0.0551,0.0543,0.0536,0.0529,
                                                                 0.1882,0.1779,0.1683,0.1594,0.1512,0.1436,0.1365,0.1301,0.1241,0.1186,0.1136,0.1090,0.1048,0.1010,0.0975,0.0944,0.0916,0.0890,0.0867,0.0847,0.0828,0.0811,0.0796,0.0782,0.0769,0.0757,0.0745,0.0734,0.0724,0.0714,0.0703,0.0693,0.0683,0.0673,0.0663,0.0653,0.0642,0.0632,0.0621,0.0611,0.0600,0.0590,0.0580,0.0570,0.0561,0.0552,0.0544,0.0536,0.0529,0.0522,
                                                                 0.1861,0.1757,0.1660,0.1570,0.1487,0.1411,0.1341,0.1276,0.1217,0.1162,0.1113,0.1068,0.1027,0.0990,0.0956,0.0926,0.0899,0.0875,0.0853,0.0833,0.0816,0.0799,0.0785,0.0771,0.0759,0.0747,0.0736,0.0725,0.0715,0.0705,0.0695,0.0685,0.0675,0.0665,0.0655,0.0644,0.0634,0.0623,0.0613,0.0602,0.0592,0.0582,0.0572,0.0563,0.0553,0.0545,0.0537,0.0529,0.0522,0.0516,
                                                                 0.1839,0.1734,0.1636,0.1546,0.1463,0.1386,0.1316,0.1252,0.1193,0.1139,0.1091,0.1046,0.1007,0.0971,0.0939,0.0910,0.0884,0.0860,0.0840,0.0821,0.0804,0.0788,0.0774,0.0761,0.0749,0.0738,0.0727,0.0717,0.0707,0.0697,0.0687,0.0677,0.0667,0.0657,0.0646,0.0636,0.0626,0.0615,0.0605,0.0594,0.0584,0.0574,0.0564,0.0555,0.0546,0.0538,0.0530,0.0523,0.0516,0.0510,
                                                                 0.1816,0.1711,0.1612,0.1522,0.1438,0.1362,0.1292,0.1228,0.1170,0.1117,0.1069,0.1026,0.0988,0.0953,0.0922,0.0894,0.0869,0.0847,0.0827,0.0809,0.0793,0.0778,0.0765,0.0752,0.0741,0.0730,0.0719,0.0709,0.0699,0.0689,0.0679,0.0669,0.0659,0.0649,0.0639,0.0628,0.0618,0.0607,0.0597,0.0586,0.0576,0.0567,0.0557,0.0548,0.0539,0.0531,0.0524,0.0517,0.0510,0.0504,
                                                                 0.1793,0.1687,0.1588,0.1497,0.1414,0.1338,0.1268,0.1205,0.1147,0.1096,0.1049,0.1007,0.0970,0.0936,0.0906,0.0879,0.0856,0.0834,0.0815,0.0798,0.0783,0.0769,0.0756,0.0743,0.0732,0.0721,0.0711,0.0701,0.0691,0.0681,0.0671,0.0661,0.0651,0.0641,0.0631,0.0620,0.0610,0.0599,0.0589,0.0579,0.0569,0.0559,0.0550,0.0541,0.0533,0.0525,0.0518,0.0511,0.0505,0.0499,
                                                                 0.1769,0.1663,0.1564,0.1473,0.1390,0.1314,0.1245,0.1182,0.1126,0.1075,0.1030,0.0989,0.0953,0.0920,0.0891,0.0866,0.0843,0.0823,0.0804,0.0788,0.0773,0.0759,0.0747,0.0735,0.0724,0.0714,0.0703,0.0693,0.0684,0.0674,0.0664,0.0654,0.0644,0.0634,0.0623,0.0613,0.0602,0.0592,0.0582,0.0572,0.0562,0.0552,0.0543,0.0535,0.0527,0.0519,0.0512,0.0506,0.0500,0.0494,
                                                                 0.1745,0.1638,0.1539,0.1448,0.1366,0.1290,0.1222,0.1160,0.1105,0.1056,0.1011,0.0972,0.0937,0.0905,0.0878,0.0853,0.0831,0.0812,0.0794,0.0778,0.0764,0.0751,0.0739,0.0727,0.0716,0.0706,0.0696,0.0686,0.0676,0.0666,0.0657,0.0647,0.0636,0.0626,0.0616,0.0606,0.0595,0.0585,0.0575,0.0565,0.0555,0.0546,0.0537,0.0529,0.0521,0.0514,0.0507,0.0501,0.0495,0.0490,
                                                                 0.1721,0.1613,0.1515,0.1424,0.1342,0.1267,0.1200,0.1140,0.1085,0.1037,0.0994,0.0956,0.0922,0.0891,0.0865,0.0841,0.0820,0.0801,0.0784,0.0769,0.0755,0.0742,0.0731,0.0720,0.0709,0.0699,0.0689,0.0679,0.0669,0.0659,0.0650,0.0640,0.0629,0.0619,0.0609,0.0598,0.0588,0.0578,0.0568,0.0558,0.0549,0.0540,0.0531,0.0523,0.0516,0.0509,0.0502,0.0496,0.0491,0.0486,
                                                                 0.1696,0.1589,0.1490,0.1400,0.1319,0.1245,0.1179,0.1119,0.1066,0.1019,0.0977,0.0940,0.0907,0.0878,0.0853,0.0830,0.0809,0.0791,0.0775,0.0760,0.0747,0.0735,0.0723,0.0712,0.0702,0.0692,0.0682,0.0672,0.0662,0.0653,0.0643,0.0633,0.0623,0.0612,0.0602,0.0592,0.0582,0.0571,0.0562,0.0552,0.0543,0.0534,0.0526,0.0518,0.0511,0.0504,0.0498,0.0492,0.0487,0.0483,
                                                                 0.1671,0.1564,0.1466,0.1377,0.1296,0.1224,0.1159,0.1100,0.1049,0.1003,0.0962,0.0926,0.0894,0.0866,0.0841,0.0819,0.0800,0.0782,0.0767,0.0752,0.0739,0.0727,0.0716,0.0705,0.0695,0.0685,0.0675,0.0666,0.0656,0.0646,0.0636,0.0626,0.0616,0.0606,0.0596,0.0585,0.0575,0.0565,0.0556,0.0546,0.0537,0.0529,0.0521,0.0513,0.0506,0.0500,0.0494,0.0489,0.0484,0.0480,
                                                                 0.1646,0.1540,0.1443,0.1354,0.1275,0.1203,0.1139,0.1082,0.1032,0.0987,0.0947,0.0912,0.0882,0.0855,0.0831,0.0809,0.0791,0.0774,0.0759,0.0745,0.0732,0.0720,0.0709,0.0699,0.0689,0.0679,0.0669,0.0659,0.0650,0.0640,0.0630,0.0620,0.0610,0.0599,0.0589,0.0579,0.0569,0.0559,0.0550,0.0541,0.0532,0.0524,0.0516,0.0509,0.0502,0.0496,0.0490,0.0485,0.0481,0.0477,
                                                                 0.1622,0.1516,0.1420,0.1332,0.1254,0.1184,0.1121,0.1065,0.1016,0.0972,0.0934,0.0900,0.0870,0.0844,0.0821,0.0800,0.0782,0.0766,0.0751,0.0738,0.0725,0.0714,0.0703,0.0693,0.0683,0.0673,0.0663,0.0653,0.0644,0.0634,0.0624,0.0614,0.0604,0.0593,0.0583,0.0573,0.0563,0.0554,0.0544,0.0536,0.0527,0.0519,0.0512,0.0505,0.0498,0.0492,0.0487,0.0482,0.0478,0.0474,
                                                                 0.1598,0.1493,0.1397,0.1311,0.1234,0.1165,0.1103,0.1049,0.1001,0.0959,0.0921,0.0888,0.0860,0.0834,0.0812,0.0792,0.0774,0.0758,0.0744,0.0731,0.0719,0.0708,0.0697,0.0687,0.0677,0.0667,0.0657,0.0648,0.0638,0.0628,0.0618,0.0608,0.0598,0.0588,0.0578,0.0568,0.0558,0.0548,0.0539,0.0531,0.0522,0.0515,0.0508,0.0501,0.0495,0.0489,0.0484,0.0480,0.0475,0.0472,
                                                                 0.1574,0.1470,0.1376,0.1291,0.1215,0.1147,0.1087,0.1034,0.0987,0.0946,0.0910,0.0878,0.0850,0.0825,0.0803,0.0784,0.0767,0.0752,0.0738,0.0725,0.0713,0.0702,0.0691,0.0681,0.0671,0.0662,0.0652,0.0642,0.0632,0.0622,0.0612,0.0602,0.0592,0.0582,0.0572,0.0562,0.0553,0.0543,0.0535,0.0526,0.0518,0.0511,0.0504,0.0497,0.0492,0.0486,0.0481,0.0477,0.0473,0.0470,
                                                                 0.1552,0.1449,0.1356,0.1272,0.1197,0.1131,0.1072,0.1020,0.0974,0.0934,0.0899,0.0868,0.0841,0.0817,0.0796,0.0777,0.0760,0.0745,0.0732,0.0719,0.0708,0.0697,0.0686,0.0676,0.0666,0.0656,0.0647,0.0637,0.0627,0.0617,0.0607,0.0597,0.0587,0.0577,0.0567,0.0557,0.0548,0.0539,0.0530,0.0522,0.0514,0.0507,0.0500,0.0494,0.0488,0.0483,0.0479,0.0475,0.0471,0.0468,
                                                                 0.1530,0.1428,0.1337,0.1254,0.1181,0.1116,0.1058,0.1007,0.0963,0.0923,0.0889,0.0859,0.0833,0.0809,0.0789,0.0770,0.0754,0.0740,0.0726,0.0714,0.0702,0.0692,0.0681,0.0671,0.0661,0.0651,0.0642,0.0632,0.0622,0.0612,0.0602,0.0592,0.0582,0.0572,0.0562,0.0552,0.0543,0.0534,0.0526,0.0518,0.0510,0.0503,0.0497,0.0491,0.0486,0.0481,0.0476,0.0473,0.0469,0.0466,
                                                                 0.1509,0.1409,0.1318,0.1237,0.1165,0.1102,0.1045,0.0995,0.0952,0.0914,0.0880,0.0851,0.0825,0.0802,0.0782,0.0764,0.0749,0.0734,0.0721,0.0709,0.0698,0.0687,0.0676,0.0666,0.0656,0.0647,0.0637,0.0627,0.0617,0.0607,0.0597,0.0587,0.0577,0.0567,0.0557,0.0548,0.0539,0.0530,0.0522,0.0514,0.0507,0.0500,0.0494,0.0488,0.0483,0.0478,0.0474,0.0471,0.0467,0.0464,
                                                                 0.1489,0.1391,0.1302,0.1222,0.1151,0.1088,0.1033,0.0985,0.0942,0.0905,0.0872,0.0843,0.0818,0.0796,0.0776,0.0759,0.0743,0.0729,0.0716,0.0704,0.0693,0.0682,0.0672,0.0662,0.0652,0.0642,0.0632,0.0622,0.0612,0.0602,0.0592,0.0582,0.0572,0.0562,0.0553,0.0543,0.0534,0.0526,0.0518,0.0510,0.0503,0.0497,0.0491,0.0486,0.0481,0.0476,0.0472,0.0469,0.0466,0.0463,
                                                                 0.1471,0.1374,0.1286,0.1208,0.1138,0.1077,0.1022,0.0975,0.0933,0.0897,0.0865,0.0837,0.0812,0.0790,0.0771,0.0754,0.0739,0.0725,0.0712,0.0700,0.0689,0.0678,0.0668,0.0658,0.0648,0.0638,0.0628,0.0618,0.0607,0.0597,0.0587,0.0577,0.0567,0.0558,0.0548,0.0539,0.0530,0.0522,0.0514,0.0507,0.0500,0.0494,0.0488,0.0483,0.0478,0.0474,0.0470,0.0467,0.0464,0.0462,
                                                                 0.1454,0.1358,0.1271,0.1194,0.1126,0.1066,0.1013,0.0966,0.0925,0.0889,0.0858,0.0830,0.0806,0.0785,0.0766,0.0749,0.0734,0.0720,0.0708,0.0696,0.0685,0.0674,0.0664,0.0653,0.0643,0.0633,0.0623,0.0613,0.0603,0.0593,0.0583,0.0573,0.0563,0.0553,0.0544,0.0535,0.0526,0.0518,0.0511,0.0504,0.0497,0.0491,0.0486,0.0481,0.0476,0.0472,0.0469,0.0466,0.0463,0.0460,
                                                                 0.1438,0.1343,0.1258,0.1182,0.1115,0.1056,0.1004,0.0958,0.0918,0.0883,0.0852,0.0825,0.0801,0.0780,0.0762,0.0745,0.0730,0.0716,0.0704,0.0692,0.0681,0.0670,0.0660,0.0649,0.0639,0.0629,0.0619,0.0609,0.0598,0.0588,0.0578,0.0568,0.0558,0.0549,0.0540,0.0531,0.0523,0.0515,0.0507,0.0500,0.0494,0.0488,0.0483,0.0479,0.0474,0.0471,0.0467,0.0464,0.0462,0.0459,
                                                                 0.1423,0.1330,0.1246,0.1171,0.1105,0.1047,0.0996,0.0951,0.0911,0.0877,0.0846,0.0820,0.0797,0.0776,0.0758,0.0741,0.0726,0.0713,0.0700,0.0688,0.0677,0.0666,0.0656,0.0645,0.0635,0.0625,0.0615,0.0604,0.0594,0.0584,0.0574,0.0564,0.0554,0.0545,0.0536,0.0527,0.0519,0.0511,0.0504,0.0498,0.0492,0.0486,0.0481,0.0477,0.0472,0.0469,0.0466,0.0463,0.0460,0.0458,
                                                                 0.1410,0.1317,0.1235,0.1161,0.1096,0.1039,0.0988,0.0944,0.0905,0.0871,0.0841,0.0815,0.0792,0.0772,0.0754,0.0737,0.0723,0.0709,0.0697,0.0685,0.0673,0.0663,0.0652,0.0641,0.0631,0.0621,0.0610,0.0600,0.0590,0.0580,0.0569,0.0560,0.0550,0.0541,0.0532,0.0523,0.0515,0.0508,0.0501,0.0495,0.0489,0.0484,0.0479,0.0475,0.0471,0.0467,0.0464,0.0462,0.0459,0.0457,
                                                                 0.1397,0.1306,0.1225,0.1152,0.1088,0.1031,0.0982,0.0938,0.0900,0.0866,0.0837,0.0811,0.0788,0.0768,0.0750,0.0734,0.0719,0.0706,0.0693,0.0681,0.0670,0.0659,0.0648,0.0638,0.0627,0.0617,0.0606,0.0596,0.0585,0.0575,0.0565,0.0555,0.0546,0.0537,0.0528,0.0520,0.0512,0.0505,0.0498,0.0492,0.0487,0.0481,0.0477,0.0473,0.0469,0.0466,0.0463,0.0460,0.0458,0.0456,
                                                                 0.1386,0.1296,0.1215,0.1144,0.1080,0.1025,0.0976,0.0932,0.0895,0.0862,0.0833,0.0807,0.0785,0.0765,0.0747,0.0731,0.0716,0.0703,0.0690,0.0678,0.0666,0.0655,0.0645,0.0634,0.0623,0.0613,0.0602,0.0592,0.0581,0.0571,0.0561,0.0551,0.0542,0.0533,0.0524,0.0516,0.0509,0.0502,0.0495,0.0490,0.0484,0.0479,0.0475,0.0471,0.0468,0.0465,0.0462,0.0459,0.0457,0.0455,
                                                                 0.1375,0.1286,0.1207,0.1136,0.1074,0.1019,0.0970,0.0928,0.0890,0.0858,0.0829,0.0804,0.0781,0.0761,0.0744,0.0728,0.0713,0.0699,0.0687,0.0675,0.0663,0.0652,0.0641,0.0630,0.0619,0.0609,0.0598,0.0587,0.0577,0.0567,0.0557,0.0547,0.0538,0.0529,0.0521,0.0513,0.0506,0.0499,0.0493,0.0487,0.0482,0.0477,0.0473,0.0470,0.0466,0.0463,0.0461,0.0458,0.0456,0.0455,
                                                                 0.1366,0.1278,0.1199,0.1129,0.1068,0.1013,0.0965,0.0923,0.0886,0.0854,0.0825,0.0800,0.0778,0.0758,0.0741,0.0725,0.0710,0.0696,0.0683,0.0671,0.0660,0.0648,0.0637,0.0626,0.0615,0.0605,0.0594,0.0583,0.0573,0.0563,0.0553,0.0543,0.0534,0.0526,0.0518,0.0510,0.0503,0.0496,0.0490,0.0485,0.0480,0.0475,0.0472,0.0468,0.0465,0.0462,0.0460,0.0457,0.0456,0.0454,
                                                                 0.1357,0.1270,0.1192,0.1123,0.1062,0.1008,0.0961,0.0919,0.0882,0.0850,0.0822,0.0797,0.0775,0.0755,0.0738,0.0722,0.0707,0.0693,0.0680,0.0668,0.0656,0.0645,0.0633,0.0622,0.0611,0.0600,0.0590,0.0579,0.0569,0.0559,0.0549,0.0540,0.0531,0.0522,0.0514,0.0507,0.0500,0.0494,0.0488,0.0483,0.0478,0.0474,0.0470,0.0467,0.0464,0.0461,0.0459,0.0457,0.0455,0.0453,
                                                                 0.1349,0.1263,0.1186,0.1118,0.1057,0.1003,0.0956,0.0915,0.0879,0.0847,0.0819,0.0794,0.0772,0.0752,0.0735,0.0719,0.0704,0.0690,0.0677,0.0665,0.0653,0.0641,0.0630,0.0618,0.0607,0.0596,0.0586,0.0575,0.0565,0.0555,0.0545,0.0536,0.0527,0.0519,0.0511,0.0504,0.0497,0.0491,0.0486,0.0481,0.0476,0.0472,0.0468,0.0465,0.0462,0.0460,0.0458,0.0456,0.0454,0.0453,
                                                                 0.1342,0.1257,0.1180,0.1112,0.1052,0.0999,0.0953,0.0912,0.0875,0.0844,0.0816,0.0791,0.0769,0.0750,0.0732,0.0716,0.0701,0.0687,0.0674,0.0661,0.0649,0.0637,0.0626,0.0615,0.0603,0.0592,0.0582,0.0571,0.0561,0.0551,0.0541,0.0532,0.0524,0.0515,0.0508,0.0501,0.0495,0.0489,0.0483,0.0479,0.0474,0.0470,0.0467,0.0464,0.0461,0.0459,0.0457,0.0455,0.0453,0.0452,
                                                                 0.1335,0.1251,0.1175,0.1108,0.1048,0.0995,0.0949,0.0908,0.0872,0.0841,0.0813,0.0789,0.0767,0.0747,0.0729,0.0713,0.0698,0.0684,0.0671,0.0658,0.0646,0.0634,0.0622,0.0611,0.0599,0.0588,0.0578,0.0567,0.0557,0.0547,0.0538,0.0529,0.0520,0.0512,0.0505,0.0498,0.0492,0.0486,0.0481,0.0477,0.0473,0.0469,0.0466,0.0463,0.0460,0.0458,0.0456,0.0454,0.0453,0.0451};
        
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
    
    case 5:
    {
        static const std::vector<double> arrayCurvatureMin   = {0.0025,-0.0026,-0.0082,-0.0140,-0.0194,-0.0238,-0.0270,-0.0290,-0.0301,-0.0307,-0.0309,-0.0310,-0.0310,-0.0309,-0.0308,-0.0306,-0.0303,-0.0300,-0.0296,-0.0292,-0.0289,-0.0288,-0.0288,-0.0289,-0.0290,-0.0293,-0.0295,-0.0298,-0.0301,-0.0303,-0.0306,-0.0308,-0.0310,-0.0312,-0.0313,-0.0315,-0.0316,-0.0317,-0.0318,-0.0319,-0.0319,-0.0320,-0.0320,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,
                                                                0.0032,-0.0012,-0.0062,-0.0115,-0.0169,-0.0216,-0.0253,-0.0279,-0.0296,-0.0305,-0.0310,-0.0313,-0.0314,-0.0315,-0.0315,-0.0313,-0.0311,-0.0307,-0.0303,-0.0299,-0.0296,-0.0294,-0.0293,-0.0293,-0.0294,-0.0296,-0.0298,-0.0300,-0.0303,-0.0305,-0.0307,-0.0309,-0.0311,-0.0313,-0.0314,-0.0315,-0.0316,-0.0317,-0.0318,-0.0319,-0.0320,-0.0320,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,
                                                                0.0035,-0.0002,-0.0044,-0.0091,-0.0142,-0.0191,-0.0233,-0.0264,-0.0287,-0.0301,-0.0309,-0.0315,-0.0318,-0.0320,-0.0321,-0.0320,-0.0318,-0.0315,-0.0310,-0.0306,-0.0302,-0.0299,-0.0297,-0.0297,-0.0298,-0.0299,-0.0301,-0.0303,-0.0305,-0.0307,-0.0309,-0.0310,-0.0312,-0.0313,-0.0315,-0.0316,-0.0317,-0.0318,-0.0319,-0.0319,-0.0320,-0.0321,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,
                                                                0.0035,0.0005 ,-0.0029,-0.0069,-0.0115,-0.0163,-0.0208,-0.0245,-0.0274,-0.0293,-0.0307,-0.0315,-0.0321,-0.0324,-0.0326,-0.0326,-0.0325,-0.0322,-0.0317,-0.0312,-0.0308,-0.0304,-0.0302,-0.0301,-0.0301,-0.0302,-0.0303,-0.0305,-0.0306,-0.0308,-0.0310,-0.0311,-0.0313,-0.0314,-0.0315,-0.0317,-0.0318,-0.0318,-0.0319,-0.0320,-0.0320,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,
                                                                0.0034,0.0009 ,-0.0018,-0.0050,-0.0089,-0.0134,-0.0180,-0.0222,-0.0257,-0.0283,-0.0301,-0.0314,-0.0322,-0.0328,-0.0331,-0.0332,-0.0331,-0.0328,-0.0324,-0.0319,-0.0314,-0.0310,-0.0307,-0.0305,-0.0304,-0.0304,-0.0305,-0.0307,-0.0308,-0.0310,-0.0311,-0.0313,-0.0314,-0.0315,-0.0316,-0.0317,-0.0318,-0.0319,-0.0319,-0.0320,-0.0321,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,
                                                                0.0030,0.0010 ,-0.0010,-0.0035,-0.0066,-0.0105,-0.0150,-0.0195,-0.0236,-0.0268,-0.0293,-0.0310,-0.0322,-0.0330,-0.0335,-0.0337,-0.0337,-0.0334,-0.0330,-0.0325,-0.0319,-0.0315,-0.0311,-0.0309,-0.0307,-0.0307,-0.0308,-0.0308,-0.0310,-0.0311,-0.0312,-0.0314,-0.0315,-0.0316,-0.0317,-0.0318,-0.0319,-0.0319,-0.0320,-0.0320,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,
                                                                0.0026,0.0010 ,-0.0005,-0.0023,-0.0047,-0.0079,-0.0120,-0.0165,-0.0210,-0.0250,-0.0281,-0.0304,-0.0320,-0.0331,-0.0338,-0.0341,-0.0342,-0.0340,-0.0336,-0.0330,-0.0325,-0.0319,-0.0315,-0.0312,-0.0310,-0.0310,-0.0310,-0.0310,-0.0311,-0.0312,-0.0313,-0.0314,-0.0316,-0.0317,-0.0317,-0.0318,-0.0319,-0.0320,-0.0320,-0.0321,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,
                                                                0.0021,0.0009 ,-0.0002,-0.0014,-0.0031,-0.0056,-0.0092,-0.0135,-0.0182,-0.0227,-0.0265,-0.0295,-0.0316,-0.0330,-0.0340,-0.0345,-0.0346,-0.0345,-0.0341,-0.0336,-0.0330,-0.0324,-0.0319,-0.0316,-0.0313,-0.0312,-0.0312,-0.0312,-0.0313,-0.0313,-0.0314,-0.0315,-0.0316,-0.0317,-0.0318,-0.0319,-0.0319,-0.0320,-0.0320,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,
                                                                0.0016,0.0006 ,-0.0000,-0.0008,-0.0019,-0.0037,-0.0066,-0.0105,-0.0152,-0.0201,-0.0245,-0.0281,-0.0308,-0.0327,-0.0340,-0.0347,-0.0350,-0.0349,-0.0346,-0.0341,-0.0335,-0.0328,-0.0323,-0.0319,-0.0316,-0.0314,-0.0314,-0.0314,-0.0314,-0.0315,-0.0315,-0.0316,-0.0317,-0.0318,-0.0318,-0.0319,-0.0320,-0.0320,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,
                                                                0.0012,0.0004 ,0.0000 ,-0.0003,-0.0009,-0.0022,-0.0045,-0.0078,-0.0122,-0.0172,-0.0221,-0.0263,-0.0297,-0.0321,-0.0338,-0.0348,-0.0352,-0.0353,-0.0350,-0.0345,-0.0339,-0.0333,-0.0327,-0.0322,-0.0319,-0.0317,-0.0316,-0.0315,-0.0315,-0.0316,-0.0316,-0.0317,-0.0318,-0.0318,-0.0319,-0.0320,-0.0320,-0.0321,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,
                                                                0.0008,0.0001 ,0.0000 ,-0.0000,-0.0003,-0.0010,-0.0027,-0.0054,-0.0094,-0.0142,-0.0193,-0.0241,-0.0282,-0.0312,-0.0333,-0.0347,-0.0354,-0.0356,-0.0354,-0.0350,-0.0344,-0.0337,-0.0331,-0.0326,-0.0322,-0.0319,-0.0317,-0.0317,-0.0317,-0.0317,-0.0317,-0.0318,-0.0318,-0.0319,-0.0319,-0.0320,-0.0320,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,
                                                                0.0005,-0.0001,-0.0000,0.0002 ,0.0002 ,-0.0002,-0.0013,-0.0035,-0.0068,-0.0112,-0.0164,-0.0216,-0.0262,-0.0299,-0.0326,-0.0343,-0.0353,-0.0357,-0.0357,-0.0353,-0.0347,-0.0341,-0.0334,-0.0329,-0.0324,-0.0321,-0.0319,-0.0318,-0.0318,-0.0318,-0.0318,-0.0318,-0.0319,-0.0319,-0.0320,-0.0320,-0.0321,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,
                                                                0.0002,-0.0003,-0.0001,0.0003 ,0.0006 ,0.0005 ,-0.0002,-0.0019,-0.0046,-0.0085,-0.0134,-0.0187,-0.0238,-0.0282,-0.0315,-0.0337,-0.0351,-0.0358,-0.0359,-0.0356,-0.0351,-0.0344,-0.0338,-0.0332,-0.0327,-0.0323,-0.0321,-0.0320,-0.0319,-0.0319,-0.0319,-0.0319,-0.0319,-0.0320,-0.0320,-0.0321,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,
                                                                0.0001,-0.0004,-0.0001,0.0004 ,0.0008 ,0.0009 ,0.0005 ,-0.0007,-0.0028,-0.0061,-0.0105,-0.0157,-0.0211,-0.0260,-0.0299,-0.0328,-0.0346,-0.0357,-0.0360,-0.0359,-0.0354,-0.0348,-0.0341,-0.0335,-0.0329,-0.0325,-0.0323,-0.0321,-0.0320,-0.0320,-0.0319,-0.0320,-0.0320,-0.0320,-0.0321,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,
                                                                0.0002,-0.0004,-0.0002,0.0004 ,0.0010 ,0.0012 ,0.0011 ,0.0002 ,-0.0014,-0.0041,-0.0080,-0.0128,-0.0182,-0.0234,-0.0280,-0.0315,-0.0339,-0.0353,-0.0359,-0.0360,-0.0356,-0.0351,-0.0344,-0.0337,-0.0332,-0.0327,-0.0324,-0.0322,-0.0321,-0.0320,-0.0320,-0.0320,-0.0320,-0.0321,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,
                                                                0.0004,-0.0004,-0.0002,0.0005 ,0.0011 ,0.0014 ,0.0014 ,0.0009 ,-0.0004,-0.0025,-0.0058,-0.0100,-0.0152,-0.0206,-0.0256,-0.0298,-0.0328,-0.0347,-0.0357,-0.0360,-0.0358,-0.0353,-0.0346,-0.0340,-0.0334,-0.0329,-0.0326,-0.0324,-0.0322,-0.0321,-0.0321,-0.0321,-0.0321,-0.0321,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,
                                                                0.0007,-0.0003,-0.0001,0.0005 ,0.0011 ,0.0015 ,0.0016 ,0.0013 ,0.0004 ,-0.0013,-0.0039,-0.0076,-0.0123,-0.0176,-0.0229,-0.0276,-0.0313,-0.0338,-0.0352,-0.0358,-0.0358,-0.0355,-0.0349,-0.0342,-0.0336,-0.0331,-0.0328,-0.0325,-0.0323,-0.0322,-0.0322,-0.0321,-0.0321,-0.0321,-0.0322,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,
                                                                0.0012,-0.0001,-0.0000,0.0005 ,0.0011 ,0.0016 ,0.0018 ,0.0016 ,0.0009 ,-0.0004,-0.0025,-0.0056,-0.0096,-0.0146,-0.0199,-0.0250,-0.0293,-0.0325,-0.0345,-0.0355,-0.0358,-0.0355,-0.0350,-0.0344,-0.0338,-0.0333,-0.0329,-0.0326,-0.0324,-0.0323,-0.0322,-0.0322,-0.0322,-0.0322,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,
                                                                0.0020,0.0003 ,0.0001 ,0.0006 ,0.0012 ,0.0016 ,0.0018 ,0.0017 ,0.0012 ,0.0002 ,-0.0014,-0.0039,-0.0073,-0.0117,-0.0168,-0.0221,-0.0269,-0.0307,-0.0334,-0.0349,-0.0355,-0.0355,-0.0351,-0.0346,-0.0340,-0.0335,-0.0330,-0.0327,-0.0325,-0.0324,-0.0323,-0.0322,-0.0322,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0029,0.0007 ,0.0003 ,0.0006 ,0.0012 ,0.0016 ,0.0018 ,0.0018 ,0.0014 ,0.0006 ,-0.0007,-0.0026,-0.0054,-0.0092,-0.0139,-0.0191,-0.0242,-0.0286,-0.0319,-0.0340,-0.0350,-0.0353,-0.0351,-0.0347,-0.0341,-0.0336,-0.0332,-0.0328,-0.0326,-0.0324,-0.0323,-0.0323,-0.0322,-0.0322,-0.0322,-0.0322,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0041,0.0013 ,0.0006 ,0.0007 ,0.0012 ,0.0015 ,0.0018 ,0.0017 ,0.0015 ,0.0009 ,-0.0001,-0.0017,-0.0039,-0.0070,-0.0111,-0.0160,-0.0212,-0.0260,-0.0299,-0.0327,-0.0343,-0.0350,-0.0350,-0.0347,-0.0342,-0.0337,-0.0333,-0.0329,-0.0327,-0.0325,-0.0324,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0056,0.0020 ,0.0009 ,0.0008 ,0.0012 ,0.0015 ,0.0017 ,0.0017 ,0.0015 ,0.0010 ,0.0002 ,-0.0010,-0.0027,-0.0053,-0.0087,-0.0131,-0.0181,-0.0231,-0.0276,-0.0310,-0.0332,-0.0344,-0.0348,-0.0346,-0.0343,-0.0338,-0.0334,-0.0330,-0.0328,-0.0326,-0.0325,-0.0324,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0073,0.0029 ,0.0013 ,0.0010 ,0.0012 ,0.0014 ,0.0016 ,0.0016 ,0.0014 ,0.0011 ,0.0004 ,-0.0005,-0.0019,-0.0039,-0.0067,-0.0105,-0.0150,-0.0200,-0.0249,-0.0289,-0.0318,-0.0335,-0.0343,-0.0344,-0.0342,-0.0338,-0.0335,-0.0331,-0.0328,-0.0326,-0.0325,-0.0324,-0.0324,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0094,0.0040 ,0.0018 ,0.0012 ,0.0012 ,0.0014 ,0.0015 ,0.0015 ,0.0013 ,0.0010 ,0.0005 ,-0.0002,-0.0013,-0.0028,-0.0051,-0.0082,-0.0122,-0.0169,-0.0219,-0.0264,-0.0299,-0.0323,-0.0336,-0.0341,-0.0341,-0.0338,-0.0335,-0.0332,-0.0329,-0.0327,-0.0325,-0.0324,-0.0324,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0117,0.0052 ,0.0023 ,0.0014 ,0.0012 ,0.0013 ,0.0014 ,0.0014 ,0.0012 ,0.0010 ,0.0006 ,-0.0000,-0.0009,-0.0021,-0.0038,-0.0063,-0.0097,-0.0139,-0.0187,-0.0235,-0.0276,-0.0306,-0.0325,-0.0335,-0.0338,-0.0337,-0.0334,-0.0332,-0.0329,-0.0327,-0.0326,-0.0325,-0.0324,-0.0324,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0144,0.0067 ,0.0031 ,0.0017 ,0.0013 ,0.0012 ,0.0013 ,0.0012 ,0.0011 ,0.0009 ,0.0006 ,0.0001 ,-0.0006,-0.0015,-0.0029,-0.0049,-0.0076,-0.0112,-0.0157,-0.0204,-0.0249,-0.0286,-0.0311,-0.0326,-0.0333,-0.0334,-0.0333,-0.0331,-0.0329,-0.0328,-0.0326,-0.0325,-0.0324,-0.0324,-0.0324,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0173,0.0085 ,0.0039 ,0.0020 ,0.0014 ,0.0012 ,0.0012 ,0.0011 ,0.0010 ,0.0008 ,0.0005 ,0.0001 ,-0.0004,-0.0012,-0.0022,-0.0037,-0.0059,-0.0089,-0.0128,-0.0173,-0.0219,-0.0261,-0.0293,-0.0314,-0.0325,-0.0330,-0.0331,-0.0331,-0.0329,-0.0328,-0.0326,-0.0325,-0.0325,-0.0324,-0.0324,-0.0324,-0.0324,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0205,0.0105 ,0.0049 ,0.0024 ,0.0015 ,0.0012 ,0.0011 ,0.0010 ,0.0008 ,0.0007 ,0.0004 ,0.0001 ,-0.0003,-0.0009,-0.0017,-0.0029,-0.0046,-0.0070,-0.0102,-0.0143,-0.0188,-0.0233,-0.0270,-0.0298,-0.0315,-0.0324,-0.0328,-0.0329,-0.0328,-0.0327,-0.0326,-0.0325,-0.0325,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0239,0.0128 ,0.0061 ,0.0029 ,0.0016 ,0.0012 ,0.0010 ,0.0008 ,0.0007 ,0.0005 ,0.0003 ,0.0000 ,-0.0003,-0.0008,-0.0014,-0.0023,-0.0035,-0.0054,-0.0080,-0.0115,-0.0157,-0.0202,-0.0244,-0.0278,-0.0301,-0.0315,-0.0323,-0.0326,-0.0327,-0.0326,-0.0326,-0.0325,-0.0325,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0274,0.0154 ,0.0075 ,0.0036 ,0.0019 ,0.0012 ,0.0009 ,0.0007 ,0.0006 ,0.0004 ,0.0002 ,-0.0000,-0.0003,-0.0007,-0.0012,-0.0018,-0.0028,-0.0042,-0.0063,-0.0091,-0.0128,-0.0171,-0.0215,-0.0254,-0.0283,-0.0303,-0.0315,-0.0321,-0.0324,-0.0325,-0.0325,-0.0325,-0.0325,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0310,0.0182 ,0.0092 ,0.0043 ,0.0021 ,0.0012 ,0.0008 ,0.0006 ,0.0004 ,0.0003 ,0.0001 ,-0.0001,-0.0004,-0.0007,-0.0010,-0.0015,-0.0022,-0.0033,-0.0049,-0.0071,-0.0102,-0.0141,-0.0184,-0.0226,-0.0261,-0.0288,-0.0305,-0.0315,-0.0320,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0346,0.0213 ,0.0111 ,0.0052 ,0.0025 ,0.0013 ,0.0008 ,0.0005 ,0.0003 ,0.0001 ,-0.0000,-0.0002,-0.0004,-0.0007,-0.0010,-0.0013,-0.0019,-0.0026,-0.0038,-0.0055,-0.0080,-0.0113,-0.0153,-0.0196,-0.0235,-0.0268,-0.0291,-0.0306,-0.0315,-0.0320,-0.0322,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0381,0.0246 ,0.0133 ,0.0063 ,0.0029 ,0.0014 ,0.0007 ,0.0004 ,0.0002 ,0.0000 ,-0.0001,-0.0003,-0.0005,-0.0007,-0.0009,-0.0012,-0.0016,-0.0022,-0.0030,-0.0043,-0.0062,-0.0089,-0.0124,-0.0165,-0.0207,-0.0244,-0.0273,-0.0294,-0.0307,-0.0315,-0.0319,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0415,0.0280 ,0.0157 ,0.0077 ,0.0035 ,0.0016 ,0.0007 ,0.0003 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0006,-0.0007,-0.0009,-0.0011,-0.0014,-0.0018,-0.0024,-0.0033,-0.0048,-0.0069,-0.0098,-0.0134,-0.0176,-0.0216,-0.0251,-0.0278,-0.0296,-0.0308,-0.0315,-0.0319,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0446,0.0316 ,0.0184 ,0.0092 ,0.0042 ,0.0018 ,0.0008 ,0.0003 ,-0.0000,-0.0002,-0.0004,-0.0005,-0.0006,-0.0008,-0.0009,-0.0011,-0.0013,-0.0016,-0.0020,-0.0027,-0.0037,-0.0052,-0.0075,-0.0107,-0.0145,-0.0186,-0.0226,-0.0258,-0.0282,-0.0299,-0.0309,-0.0315,-0.0319,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0475,0.0351 ,0.0214 ,0.0110 ,0.0050 ,0.0022 ,0.0008 ,0.0002 ,-0.0001,-0.0003,-0.0005,-0.0006,-0.0007,-0.0009,-0.0010,-0.0011,-0.0013,-0.0015,-0.0017,-0.0022,-0.0029,-0.0040,-0.0057,-0.0082,-0.0116,-0.0155,-0.0197,-0.0234,-0.0265,-0.0287,-0.0301,-0.0310,-0.0316,-0.0319,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0501,0.0385 ,0.0246 ,0.0131 ,0.0061 ,0.0026 ,0.0010 ,0.0002 ,-0.0002,-0.0004,-0.0005,-0.0007,-0.0008,-0.0009,-0.0010,-0.0011,-0.0012,-0.0014,-0.0015,-0.0018,-0.0023,-0.0030,-0.0042,-0.0062,-0.0089,-0.0125,-0.0166,-0.0207,-0.0243,-0.0271,-0.0291,-0.0304,-0.0312,-0.0317,-0.0320,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0524,0.0418 ,0.0279 ,0.0154 ,0.0073 ,0.0031 ,0.0011 ,0.0002 ,-0.0002,-0.0005,-0.0006,-0.0008,-0.0009,-0.0010,-0.0011,-0.0012,-0.0012,-0.0013,-0.0014,-0.0016,-0.0018,-0.0023,-0.0031,-0.0045,-0.0067,-0.0097,-0.0135,-0.0177,-0.0217,-0.0251,-0.0277,-0.0294,-0.0306,-0.0313,-0.0318,-0.0320,-0.0322,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0544,0.0449 ,0.0314 ,0.0180 ,0.0088 ,0.0038 ,0.0014 ,0.0003 ,-0.0003,-0.0005,-0.0007,-0.0008,-0.0010,-0.0010,-0.0011,-0.0012,-0.0013,-0.0013,-0.0014,-0.0014,-0.0016,-0.0018,-0.0023,-0.0032,-0.0048,-0.0072,-0.0105,-0.0146,-0.0188,-0.0227,-0.0259,-0.0282,-0.0298,-0.0308,-0.0314,-0.0318,-0.0321,-0.0322,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0562,0.0478 ,0.0348 ,0.0209 ,0.0105 ,0.0046 ,0.0017 ,0.0003 ,-0.0003,-0.0006,-0.0008,-0.0009,-0.0010,-0.0011,-0.0012,-0.0012,-0.0013,-0.0013,-0.0013,-0.0013,-0.0014,-0.0015,-0.0017,-0.0023,-0.0033,-0.0051,-0.0078,-0.0115,-0.0157,-0.0200,-0.0237,-0.0266,-0.0287,-0.0301,-0.0310,-0.0316,-0.0319,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0577,0.0503 ,0.0382 ,0.0240 ,0.0125 ,0.0056 ,0.0021 ,0.0005 ,-0.0003,-0.0006,-0.0008,-0.0010,-0.0011,-0.0012,-0.0012,-0.0013,-0.0013,-0.0013,-0.0013,-0.0013,-0.0013,-0.0012,-0.0013,-0.0016,-0.0022,-0.0034,-0.0055,-0.0086,-0.0126,-0.0170,-0.0211,-0.0247,-0.0273,-0.0292,-0.0304,-0.0312,-0.0317,-0.0320,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0590,0.0526 ,0.0415 ,0.0273 ,0.0148 ,0.0068 ,0.0026 ,0.0007 ,-0.0002,-0.0007,-0.0009,-0.0010,-0.0011,-0.0012,-0.0013,-0.0013,-0.0013,-0.0013,-0.0013,-0.0013,-0.0012,-0.0011,-0.0010,-0.0011,-0.0014,-0.0021,-0.0036,-0.0061,-0.0095,-0.0138,-0.0182,-0.0223,-0.0256,-0.0280,-0.0296,-0.0307,-0.0314,-0.0318,-0.0320,-0.0322,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0600,0.0546 ,0.0446 ,0.0307 ,0.0173 ,0.0082 ,0.0032 ,0.0009 ,-0.0002,-0.0007,-0.0009,-0.0011,-0.0012,-0.0013,-0.0013,-0.0014,-0.0014,-0.0014,-0.0013,-0.0013,-0.0012,-0.0010,-0.0009,-0.0008,-0.0008,-0.0012,-0.0021,-0.0039,-0.0068,-0.0106,-0.0151,-0.0195,-0.0234,-0.0264,-0.0286,-0.0300,-0.0309,-0.0315,-0.0319,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0609,0.0564 ,0.0474 ,0.0341 ,0.0202 ,0.0099 ,0.0040 ,0.0012 ,-0.0001,-0.0007,-0.0010,-0.0011,-0.0012,-0.0013,-0.0014,-0.0014,-0.0014,-0.0014,-0.0014,-0.0013,-0.0012,-0.0010,-0.0008,-0.0006,-0.0004,-0.0005,-0.0010,-0.0022,-0.0044,-0.0077,-0.0119,-0.0165,-0.0208,-0.0245,-0.0272,-0.0291,-0.0303,-0.0311,-0.0316,-0.0319,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0617,0.0578 ,0.0500 ,0.0375 ,0.0232 ,0.0118 ,0.0050 ,0.0016 ,0.0001 ,-0.0007,-0.0010,-0.0012,-0.0013,-0.0014,-0.0014,-0.0014,-0.0015,-0.0014,-0.0014,-0.0013,-0.0012,-0.0010,-0.0007,-0.0005,-0.0002,-0.0000,-0.0001,-0.0008,-0.0024,-0.0050,-0.0088,-0.0133,-0.0179,-0.0221,-0.0254,-0.0279,-0.0296,-0.0306,-0.0313,-0.0317,-0.0320,-0.0322,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0623,0.0591 ,0.0523 ,0.0408 ,0.0265 ,0.0140 ,0.0062 ,0.0021 ,0.0002 ,-0.0006,-0.0010,-0.0012,-0.0013,-0.0014,-0.0015,-0.0015,-0.0015,-0.0015,-0.0014,-0.0013,-0.0012,-0.0010,-0.0007,-0.0004,-0.0000,0.0003 ,0.0004 ,0.0002 ,-0.0008,-0.0028,-0.0059,-0.0101,-0.0148,-0.0194,-0.0233,-0.0264,-0.0285,-0.0300,-0.0309,-0.0315,-0.0318,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0628,0.0601 ,0.0543 ,0.0439 ,0.0299 ,0.0166 ,0.0076 ,0.0028 ,0.0005 ,-0.0005,-0.0010,-0.0012,-0.0014,-0.0014,-0.0015,-0.0015,-0.0015,-0.0015,-0.0015,-0.0014,-0.0012,-0.0010,-0.0008,-0.0004,0.0000 ,0.0005 ,0.0008 ,0.0009 ,0.0004 ,-0.0010,-0.0034,-0.0070,-0.0115,-0.0163,-0.0208,-0.0244,-0.0272,-0.0291,-0.0303,-0.0311,-0.0316,-0.0319,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,-0.0324,
                                                                0.0632,0.0610 ,0.0561 ,0.0468 ,0.0333 ,0.0194 ,0.0092 ,0.0035 ,0.0008 ,-0.0004,-0.0010,-0.0013,-0.0014,-0.0015,-0.0015,-0.0016,-0.0016,-0.0015,-0.0015,-0.0014,-0.0013,-0.0011,-0.0008,-0.0004,0.0000 ,0.0006 ,0.0010 ,0.0014 ,0.0012 ,0.0004 ,-0.0014,-0.0043,-0.0083,-0.0131,-0.0179,-0.0221,-0.0255,-0.0279,-0.0296,-0.0306,-0.0313,-0.0317,-0.0320,-0.0322,-0.0323,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,
                                                                0.0635,0.0617 ,0.0576 ,0.0494 ,0.0367 ,0.0224 ,0.0111 ,0.0045 ,0.0012 ,-0.0003,-0.0010,-0.0013,-0.0014,-0.0015,-0.0016,-0.0016,-0.0016,-0.0016,-0.0015,-0.0015,-0.0013,-0.0011,-0.0009,-0.0005,0.0000 ,0.0006 ,0.0012 ,0.0017 ,0.0019 ,0.0015 ,0.0003 ,-0.0020,-0.0054,-0.0098,-0.0147,-0.0194,-0.0234,-0.0264,-0.0286,-0.0300,-0.0309,-0.0315,-0.0318,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324,-0.0324,
                                                                0.0638,0.0623 ,0.0589 ,0.0518 ,0.0400 ,0.0256 ,0.0133 ,0.0056 ,0.0017 ,-0.0001,-0.0009,-0.0013,-0.0014,-0.0015,-0.0016,-0.0016,-0.0016,-0.0016,-0.0016,-0.0015,-0.0014,-0.0012,-0.0009,-0.0005,-0.0001,0.0005 ,0.0012 ,0.0018 ,0.0023 ,0.0023 ,0.0016 ,-0.0000,-0.0028,-0.0067,-0.0114,-0.0163,-0.0208,-0.0245,-0.0272,-0.0291,-0.0304,-0.0311,-0.0316,-0.0319,-0.0321,-0.0322,-0.0323,-0.0323,-0.0324,-0.0324};
        
        return arrayCurvatureMin[y * arrayWidth + x];
    }
    case 6:
    {
        static const std::vector<double> arrayCurvatureMin   = {0.0044,0.0029,0.0014,-0.0001,-0.0017,-0.0034,-0.0052,-0.0071,-0.0090,-0.0109,-0.0129,-0.0148,-0.0166,-0.0183,-0.0198,-0.0211,-0.0222,-0.0230,-0.0235,-0.0237,-0.0237,-0.0233,-0.0229,-0.0223,-0.0217,-0.0211,-0.0206,-0.0203,-0.0200,-0.0198,-0.0197,-0.0197,-0.0197,-0.0198,-0.0198,-0.0199,-0.0200,-0.0202,-0.0203,-0.0204,-0.0205,-0.0206,-0.0207,-0.0208,-0.0209,-0.0210,-0.0211,-0.0212,-0.0213,-0.0214,
                                                                0.0047,0.0032,0.0018,0.0003 ,-0.0012,-0.0028,-0.0045,-0.0063,-0.0082,-0.0101,-0.0120,-0.0140,-0.0159,-0.0176,-0.0193,-0.0207,-0.0220,-0.0229,-0.0236,-0.0240,-0.0241,-0.0240,-0.0236,-0.0231,-0.0225,-0.0219,-0.0214,-0.0210,-0.0207,-0.0204,-0.0203,-0.0203,-0.0203,-0.0203,-0.0204,-0.0204,-0.0205,-0.0206,-0.0207,-0.0209,-0.0210,-0.0211,-0.0212,-0.0213,-0.0214,-0.0214,-0.0215,-0.0216,-0.0217,-0.0218,
                                                                0.0049,0.0034,0.0020,0.0006 ,-0.0008,-0.0023,-0.0039,-0.0056,-0.0074,-0.0093,-0.0112,-0.0131,-0.0151,-0.0169,-0.0186,-0.0202,-0.0216,-0.0227,-0.0236,-0.0242,-0.0245,-0.0245,-0.0242,-0.0238,-0.0232,-0.0227,-0.0221,-0.0217,-0.0213,-0.0211,-0.0209,-0.0208,-0.0208,-0.0208,-0.0208,-0.0209,-0.0210,-0.0211,-0.0212,-0.0213,-0.0214,-0.0215,-0.0215,-0.0216,-0.0217,-0.0218,-0.0219,-0.0220,-0.0220,-0.0221,
                                                                0.0051,0.0036,0.0022,0.0009 ,-0.0005,-0.0019,-0.0034,-0.0050,-0.0067,-0.0085,-0.0104,-0.0123,-0.0142,-0.0161,-0.0179,-0.0196,-0.0211,-0.0224,-0.0234,-0.0242,-0.0246,-0.0248,-0.0247,-0.0244,-0.0239,-0.0234,-0.0228,-0.0223,-0.0219,-0.0216,-0.0214,-0.0213,-0.0212,-0.0212,-0.0213,-0.0213,-0.0214,-0.0215,-0.0215,-0.0216,-0.0217,-0.0218,-0.0219,-0.0220,-0.0220,-0.0221,-0.0222,-0.0222,-0.0223,-0.0224,
                                                                0.0052,0.0037,0.0024,0.0011 ,-0.0002,-0.0016,-0.0030,-0.0045,-0.0061,-0.0078,-0.0096,-0.0114,-0.0133,-0.0152,-0.0171,-0.0188,-0.0205,-0.0219,-0.0231,-0.0240,-0.0246,-0.0250,-0.0250,-0.0248,-0.0245,-0.0240,-0.0234,-0.0229,-0.0225,-0.0222,-0.0219,-0.0217,-0.0216,-0.0216,-0.0216,-0.0217,-0.0217,-0.0218,-0.0219,-0.0219,-0.0220,-0.0221,-0.0222,-0.0222,-0.0223,-0.0224,-0.0224,-0.0225,-0.0225,-0.0226,
                                                                0.0054,0.0038,0.0024,0.0012 ,-0.0000,-0.0013,-0.0026,-0.0040,-0.0055,-0.0071,-0.0088,-0.0106,-0.0124,-0.0143,-0.0162,-0.0180,-0.0197,-0.0212,-0.0226,-0.0237,-0.0245,-0.0250,-0.0252,-0.0252,-0.0249,-0.0245,-0.0240,-0.0235,-0.0230,-0.0226,-0.0223,-0.0221,-0.0220,-0.0220,-0.0219,-0.0220,-0.0220,-0.0221,-0.0221,-0.0222,-0.0223,-0.0223,-0.0224,-0.0225,-0.0225,-0.0226,-0.0226,-0.0227,-0.0227,-0.0228,
                                                                0.0055,0.0039,0.0025,0.0013 ,0.0001 ,-0.0010,-0.0022,-0.0035,-0.0049,-0.0064,-0.0080,-0.0098,-0.0115,-0.0134,-0.0152,-0.0171,-0.0188,-0.0205,-0.0219,-0.0232,-0.0241,-0.0248,-0.0252,-0.0253,-0.0252,-0.0249,-0.0244,-0.0239,-0.0235,-0.0230,-0.0227,-0.0225,-0.0223,-0.0222,-0.0222,-0.0222,-0.0222,-0.0223,-0.0224,-0.0224,-0.0225,-0.0225,-0.0226,-0.0227,-0.0227,-0.0228,-0.0228,-0.0228,-0.0229,-0.0229,
                                                                0.0057,0.0040,0.0026,0.0014 ,0.0003 ,-0.0008,-0.0020,-0.0032,-0.0044,-0.0058,-0.0073,-0.0089,-0.0106,-0.0124,-0.0142,-0.0161,-0.0179,-0.0196,-0.0211,-0.0225,-0.0236,-0.0245,-0.0251,-0.0254,-0.0254,-0.0251,-0.0248,-0.0243,-0.0238,-0.0234,-0.0231,-0.0228,-0.0226,-0.0225,-0.0224,-0.0224,-0.0224,-0.0225,-0.0225,-0.0226,-0.0226,-0.0227,-0.0228,-0.0228,-0.0229,-0.0229,-0.0229,-0.0230,-0.0230,-0.0230,
                                                                0.0059,0.0041,0.0027,0.0015 ,0.0004 ,-0.0006,-0.0017,-0.0028,-0.0040,-0.0053,-0.0066,-0.0081,-0.0098,-0.0115,-0.0132,-0.0150,-0.0168,-0.0186,-0.0202,-0.0217,-0.0230,-0.0240,-0.0247,-0.0252,-0.0254,-0.0253,-0.0250,-0.0246,-0.0242,-0.0237,-0.0233,-0.0230,-0.0228,-0.0227,-0.0226,-0.0226,-0.0226,-0.0226,-0.0227,-0.0227,-0.0228,-0.0228,-0.0229,-0.0229,-0.0230,-0.0230,-0.0231,-0.0231,-0.0231,-0.0231,
                                                                0.0063,0.0043,0.0028,0.0016 ,0.0005 ,-0.0005,-0.0015,-0.0025,-0.0035,-0.0047,-0.0060,-0.0074,-0.0089,-0.0105,-0.0122,-0.0139,-0.0157,-0.0175,-0.0192,-0.0207,-0.0221,-0.0233,-0.0242,-0.0249,-0.0252,-0.0253,-0.0251,-0.0248,-0.0244,-0.0240,-0.0236,-0.0233,-0.0230,-0.0229,-0.0228,-0.0227,-0.0227,-0.0227,-0.0228,-0.0228,-0.0229,-0.0229,-0.0230,-0.0230,-0.0231,-0.0231,-0.0231,-0.0232,-0.0232,-0.0232,
                                                                0.0067,0.0046,0.0030,0.0017 ,0.0007 ,-0.0003,-0.0012,-0.0022,-0.0032,-0.0042,-0.0054,-0.0067,-0.0081,-0.0096,-0.0112,-0.0128,-0.0146,-0.0163,-0.0181,-0.0197,-0.0212,-0.0225,-0.0236,-0.0244,-0.0249,-0.0251,-0.0251,-0.0249,-0.0245,-0.0241,-0.0238,-0.0234,-0.0232,-0.0230,-0.0229,-0.0228,-0.0228,-0.0228,-0.0228,-0.0229,-0.0229,-0.0230,-0.0230,-0.0231,-0.0231,-0.0232,-0.0232,-0.0232,-0.0232,-0.0233,
                                                                0.0072,0.0050,0.0032,0.0019 ,0.0008 ,-0.0001,-0.0010,-0.0019,-0.0028,-0.0038,-0.0048,-0.0060,-0.0073,-0.0087,-0.0102,-0.0118,-0.0134,-0.0152,-0.0169,-0.0186,-0.0201,-0.0215,-0.0228,-0.0237,-0.0244,-0.0248,-0.0250,-0.0249,-0.0246,-0.0243,-0.0239,-0.0236,-0.0233,-0.0231,-0.0229,-0.0229,-0.0228,-0.0228,-0.0229,-0.0229,-0.0230,-0.0230,-0.0231,-0.0231,-0.0232,-0.0232,-0.0232,-0.0233,-0.0233,-0.0233,
                                                                0.0079,0.0054,0.0035,0.0021 ,0.0010 ,0.0001 ,-0.0008,-0.0016,-0.0024,-0.0033,-0.0043,-0.0053,-0.0065,-0.0078,-0.0092,-0.0107,-0.0123,-0.0140,-0.0157,-0.0174,-0.0190,-0.0205,-0.0218,-0.0229,-0.0238,-0.0244,-0.0247,-0.0247,-0.0246,-0.0243,-0.0240,-0.0236,-0.0233,-0.0231,-0.0230,-0.0229,-0.0228,-0.0229,-0.0229,-0.0229,-0.0230,-0.0230,-0.0231,-0.0231,-0.0232,-0.0232,-0.0232,-0.0233,-0.0233,-0.0233,
                                                                0.0086,0.0059,0.0039,0.0024 ,0.0012 ,0.0003 ,-0.0006,-0.0013,-0.0021,-0.0029,-0.0038,-0.0048,-0.0058,-0.0070,-0.0083,-0.0097,-0.0112,-0.0128,-0.0144,-0.0161,-0.0178,-0.0193,-0.0208,-0.0220,-0.0230,-0.0238,-0.0242,-0.0244,-0.0244,-0.0242,-0.0240,-0.0236,-0.0234,-0.0231,-0.0230,-0.0229,-0.0228,-0.0228,-0.0229,-0.0229,-0.0230,-0.0230,-0.0231,-0.0231,-0.0232,-0.0232,-0.0232,-0.0233,-0.0233,-0.0233,
                                                                0.0095,0.0066,0.0044,0.0027 ,0.0015 ,0.0005 ,-0.0003,-0.0011,-0.0018,-0.0026,-0.0033,-0.0042,-0.0052,-0.0062,-0.0074,-0.0087,-0.0101,-0.0116,-0.0132,-0.0149,-0.0165,-0.0181,-0.0196,-0.0210,-0.0221,-0.0230,-0.0237,-0.0240,-0.0242,-0.0241,-0.0239,-0.0236,-0.0233,-0.0231,-0.0229,-0.0228,-0.0228,-0.0228,-0.0228,-0.0228,-0.0229,-0.0230,-0.0230,-0.0231,-0.0231,-0.0232,-0.0232,-0.0233,-0.0233,-0.0233,
                                                                0.0106,0.0073,0.0049,0.0031 ,0.0018 ,0.0007 ,-0.0001,-0.0008,-0.0015,-0.0022,-0.0029,-0.0037,-0.0046,-0.0055,-0.0066,-0.0078,-0.0091,-0.0105,-0.0120,-0.0136,-0.0152,-0.0168,-0.0184,-0.0198,-0.0211,-0.0222,-0.0230,-0.0235,-0.0238,-0.0238,-0.0237,-0.0235,-0.0233,-0.0230,-0.0229,-0.0228,-0.0227,-0.0227,-0.0227,-0.0228,-0.0228,-0.0229,-0.0230,-0.0230,-0.0231,-0.0231,-0.0232,-0.0232,-0.0233,-0.0233,
                                                                0.0118,0.0082,0.0055,0.0035 ,0.0021 ,0.0010 ,0.0001 ,-0.0006,-0.0013,-0.0019,-0.0025,-0.0033,-0.0040,-0.0049,-0.0059,-0.0070,-0.0082,-0.0095,-0.0109,-0.0124,-0.0140,-0.0156,-0.0171,-0.0186,-0.0200,-0.0212,-0.0221,-0.0228,-0.0233,-0.0235,-0.0235,-0.0233,-0.0231,-0.0229,-0.0228,-0.0227,-0.0226,-0.0226,-0.0226,-0.0227,-0.0227,-0.0228,-0.0229,-0.0230,-0.0230,-0.0231,-0.0231,-0.0232,-0.0232,-0.0233,
                                                                0.0131,0.0092,0.0062,0.0041 ,0.0025 ,0.0013 ,0.0004 ,-0.0003,-0.0010,-0.0016,-0.0022,-0.0028,-0.0035,-0.0043,-0.0052,-0.0062,-0.0073,-0.0085,-0.0098,-0.0112,-0.0127,-0.0143,-0.0159,-0.0174,-0.0188,-0.0201,-0.0212,-0.0220,-0.0226,-0.0230,-0.0231,-0.0231,-0.0229,-0.0228,-0.0226,-0.0225,-0.0225,-0.0224,-0.0225,-0.0225,-0.0226,-0.0227,-0.0228,-0.0229,-0.0230,-0.0230,-0.0231,-0.0231,-0.0232,-0.0232,
                                                                0.0146,0.0103,0.0071,0.0047 ,0.0029 ,0.0016 ,0.0007 ,-0.0001,-0.0007,-0.0013,-0.0019,-0.0025,-0.0031,-0.0038,-0.0046,-0.0055,-0.0064,-0.0076,-0.0088,-0.0101,-0.0115,-0.0130,-0.0146,-0.0161,-0.0176,-0.0189,-0.0202,-0.0211,-0.0219,-0.0224,-0.0226,-0.0227,-0.0227,-0.0226,-0.0224,-0.0223,-0.0223,-0.0223,-0.0223,-0.0224,-0.0225,-0.0226,-0.0227,-0.0228,-0.0229,-0.0229,-0.0230,-0.0231,-0.0231,-0.0232,
                                                                0.0162,0.0116,0.0080,0.0054 ,0.0034 ,0.0020 ,0.0010 ,0.0002 ,-0.0005,-0.0010,-0.0016,-0.0021,-0.0027,-0.0033,-0.0040,-0.0048,-0.0057,-0.0067,-0.0078,-0.0091,-0.0104,-0.0118,-0.0133,-0.0148,-0.0163,-0.0177,-0.0190,-0.0201,-0.0210,-0.0217,-0.0221,-0.0223,-0.0223,-0.0223,-0.0222,-0.0221,-0.0221,-0.0221,-0.0221,-0.0222,-0.0223,-0.0224,-0.0225,-0.0226,-0.0227,-0.0228,-0.0229,-0.0230,-0.0231,-0.0231,
                                                                0.0180,0.0130,0.0091,0.0061 ,0.0040 ,0.0024 ,0.0013 ,0.0004 ,-0.0002,-0.0008,-0.0013,-0.0018,-0.0023,-0.0029,-0.0035,-0.0042,-0.0050,-0.0059,-0.0069,-0.0081,-0.0093,-0.0107,-0.0121,-0.0135,-0.0150,-0.0165,-0.0178,-0.0190,-0.0201,-0.0208,-0.0214,-0.0217,-0.0219,-0.0219,-0.0219,-0.0218,-0.0218,-0.0218,-0.0219,-0.0220,-0.0221,-0.0222,-0.0223,-0.0225,-0.0226,-0.0227,-0.0228,-0.0229,-0.0230,-0.0231,
                                                                0.0199,0.0145,0.0103,0.0070 ,0.0046 ,0.0029 ,0.0016 ,0.0007 ,0.0000 ,-0.0006,-0.0011,-0.0015,-0.0020,-0.0025,-0.0031,-0.0037,-0.0044,-0.0052,-0.0061,-0.0072,-0.0083,-0.0096,-0.0109,-0.0123,-0.0138,-0.0152,-0.0166,-0.0179,-0.0190,-0.0199,-0.0206,-0.0211,-0.0214,-0.0215,-0.0215,-0.0215,-0.0215,-0.0215,-0.0216,-0.0217,-0.0218,-0.0220,-0.0221,-0.0223,-0.0224,-0.0226,-0.0227,-0.0228,-0.0229,-0.0230,
                                                                0.0219,0.0162,0.0116,0.0080 ,0.0053 ,0.0034 ,0.0020 ,0.0010 ,0.0003 ,-0.0003,-0.0008,-0.0013,-0.0017,-0.0022,-0.0027,-0.0032,-0.0039,-0.0046,-0.0054,-0.0063,-0.0074,-0.0085,-0.0098,-0.0111,-0.0125,-0.0139,-0.0153,-0.0167,-0.0179,-0.0189,-0.0197,-0.0203,-0.0207,-0.0210,-0.0211,-0.0211,-0.0211,-0.0212,-0.0213,-0.0214,-0.0216,-0.0217,-0.0219,-0.0221,-0.0223,-0.0224,-0.0226,-0.0227,-0.0228,-0.0229,
                                                                0.0240,0.0180,0.0130,0.0091 ,0.0062 ,0.0040 ,0.0025 ,0.0014 ,0.0005 ,-0.0001,-0.0006,-0.0010,-0.0014,-0.0019,-0.0023,-0.0028,-0.0034,-0.0040,-0.0048,-0.0056,-0.0065,-0.0076,-0.0087,-0.0100,-0.0113,-0.0127,-0.0141,-0.0154,-0.0167,-0.0178,-0.0188,-0.0195,-0.0200,-0.0204,-0.0206,-0.0207,-0.0207,-0.0208,-0.0209,-0.0211,-0.0213,-0.0215,-0.0217,-0.0219,-0.0221,-0.0222,-0.0224,-0.0226,-0.0227,-0.0228,
                                                                0.0262,0.0200,0.0146,0.0103 ,0.0071 ,0.0047 ,0.0029 ,0.0017 ,0.0008 ,0.0002 ,-0.0004,-0.0008,-0.0012,-0.0016,-0.0020,-0.0024,-0.0029,-0.0035,-0.0042,-0.0049,-0.0058,-0.0067,-0.0078,-0.0089,-0.0102,-0.0115,-0.0128,-0.0142,-0.0155,-0.0167,-0.0177,-0.0186,-0.0192,-0.0197,-0.0200,-0.0202,-0.0203,-0.0204,-0.0205,-0.0207,-0.0209,-0.0211,-0.0214,-0.0216,-0.0218,-0.0220,-0.0222,-0.0224,-0.0226,-0.0227,
                                                                0.0285,0.0220,0.0163,0.0117 ,0.0081 ,0.0054 ,0.0035 ,0.0021 ,0.0011 ,0.0004 ,-0.0001,-0.0006,-0.0010,-0.0013,-0.0017,-0.0021,-0.0026,-0.0031,-0.0036,-0.0043,-0.0051,-0.0059,-0.0069,-0.0079,-0.0091,-0.0103,-0.0116,-0.0130,-0.0143,-0.0155,-0.0166,-0.0175,-0.0183,-0.0189,-0.0193,-0.0196,-0.0198,-0.0199,-0.0201,-0.0203,-0.0205,-0.0208,-0.0210,-0.0213,-0.0216,-0.0218,-0.0220,-0.0222,-0.0224,-0.0226,
                                                                0.0308,0.0242,0.0182,0.0131 ,0.0092 ,0.0062 ,0.0041 ,0.0025 ,0.0015 ,0.0007 ,0.0001 ,-0.0004,-0.0008,-0.0011,-0.0014,-0.0018,-0.0022,-0.0027,-0.0032,-0.0037,-0.0044,-0.0052,-0.0061,-0.0070,-0.0081,-0.0092,-0.0105,-0.0118,-0.0130,-0.0143,-0.0154,-0.0165,-0.0173,-0.0180,-0.0185,-0.0189,-0.0192,-0.0194,-0.0196,-0.0198,-0.0201,-0.0204,-0.0207,-0.0210,-0.0213,-0.0215,-0.0218,-0.0220,-0.0222,-0.0224,
                                                                0.0332,0.0264,0.0201,0.0147 ,0.0104 ,0.0072 ,0.0048 ,0.0030 ,0.0018 ,0.0009 ,0.0003 ,-0.0002,-0.0005,-0.0009,-0.0012,-0.0015,-0.0019,-0.0023,-0.0027,-0.0033,-0.0039,-0.0045,-0.0053,-0.0062,-0.0072,-0.0082,-0.0094,-0.0106,-0.0118,-0.0131,-0.0143,-0.0154,-0.0163,-0.0171,-0.0177,-0.0182,-0.0185,-0.0188,-0.0191,-0.0193,-0.0196,-0.0199,-0.0203,-0.0206,-0.0209,-0.0212,-0.0215,-0.0218,-0.0220,-0.0223,
                                                                0.0355,0.0287,0.0222,0.0165 ,0.0118 ,0.0082 ,0.0055 ,0.0036 ,0.0022 ,0.0012 ,0.0006 ,0.0000 ,-0.0004,-0.0007,-0.0010,-0.0013,-0.0016,-0.0020,-0.0024,-0.0028,-0.0033,-0.0039,-0.0046,-0.0054,-0.0063,-0.0073,-0.0084,-0.0095,-0.0107,-0.0119,-0.0131,-0.0142,-0.0152,-0.0161,-0.0168,-0.0174,-0.0178,-0.0182,-0.0185,-0.0188,-0.0191,-0.0194,-0.0198,-0.0202,-0.0205,-0.0209,-0.0212,-0.0215,-0.0218,-0.0221,
                                                                0.0379,0.0310,0.0244,0.0183 ,0.0133 ,0.0093 ,0.0063 ,0.0042 ,0.0027 ,0.0016 ,0.0008 ,0.0003 ,-0.0002,-0.0005,-0.0008,-0.0011,-0.0014,-0.0017,-0.0020,-0.0024,-0.0029,-0.0034,-0.0040,-0.0047,-0.0055,-0.0064,-0.0074,-0.0085,-0.0096,-0.0108,-0.0119,-0.0131,-0.0141,-0.0151,-0.0159,-0.0165,-0.0170,-0.0175,-0.0178,-0.0182,-0.0185,-0.0189,-0.0193,-0.0197,-0.0201,-0.0205,-0.0209,-0.0212,-0.0216,-0.0218,
                                                                0.0401,0.0334,0.0266,0.0203 ,0.0149 ,0.0106 ,0.0073 ,0.0049 ,0.0031 ,0.0019 ,0.0011 ,0.0005 ,0.0000 ,-0.0003,-0.0006,-0.0009,-0.0012,-0.0014,-0.0017,-0.0021,-0.0025,-0.0030,-0.0035,-0.0041,-0.0048,-0.0056,-0.0065,-0.0075,-0.0086,-0.0097,-0.0108,-0.0119,-0.0130,-0.0140,-0.0149,-0.0156,-0.0162,-0.0167,-0.0171,-0.0175,-0.0179,-0.0183,-0.0187,-0.0192,-0.0196,-0.0201,-0.0205,-0.0209,-0.0213,-0.0216,
                                                                0.0423,0.0358,0.0289,0.0224 ,0.0167 ,0.0120 ,0.0083 ,0.0056 ,0.0037 ,0.0023 ,0.0014 ,0.0007 ,0.0002 ,-0.0002,-0.0004,-0.0007,-0.0010,-0.0012,-0.0015,-0.0018,-0.0021,-0.0025,-0.0030,-0.0036,-0.0042,-0.0049,-0.0057,-0.0066,-0.0076,-0.0086,-0.0097,-0.0108,-0.0119,-0.0129,-0.0138,-0.0146,-0.0153,-0.0159,-0.0163,-0.0168,-0.0172,-0.0177,-0.0181,-0.0186,-0.0191,-0.0196,-0.0201,-0.0205,-0.0209,-0.0213,
                                                                0.0444,0.0381,0.0313,0.0246 ,0.0185 ,0.0135 ,0.0095 ,0.0065 ,0.0043 ,0.0028 ,0.0017 ,0.0010 ,0.0004 ,0.0000 ,-0.0003,-0.0005,-0.0008,-0.0010,-0.0012,-0.0015,-0.0018,-0.0022,-0.0026,-0.0031,-0.0036,-0.0043,-0.0050,-0.0058,-0.0067,-0.0077,-0.0087,-0.0098,-0.0108,-0.0119,-0.0128,-0.0137,-0.0144,-0.0150,-0.0155,-0.0160,-0.0165,-0.0170,-0.0175,-0.0180,-0.0185,-0.0191,-0.0196,-0.0201,-0.0205,-0.0209,
                                                                0.0465,0.0404,0.0336,0.0269 ,0.0205 ,0.0151 ,0.0108 ,0.0074 ,0.0050 ,0.0033 ,0.0021 ,0.0012 ,0.0006 ,0.0002 ,-0.0001,-0.0004,-0.0006,-0.0008,-0.0010,-0.0013,-0.0015,-0.0019,-0.0022,-0.0026,-0.0031,-0.0037,-0.0043,-0.0051,-0.0059,-0.0068,-0.0077,-0.0087,-0.0098,-0.0108,-0.0118,-0.0126,-0.0134,-0.0141,-0.0147,-0.0152,-0.0157,-0.0162,-0.0168,-0.0173,-0.0179,-0.0185,-0.0191,-0.0196,-0.0201,-0.0206,
                                                                0.0483,0.0426,0.0360,0.0292 ,0.0226 ,0.0169 ,0.0121 ,0.0085 ,0.0058 ,0.0038 ,0.0025 ,0.0015 ,0.0009 ,0.0004 ,0.0000 ,-0.0002,-0.0004,-0.0006,-0.0008,-0.0011,-0.0013,-0.0016,-0.0019,-0.0022,-0.0027,-0.0032,-0.0037,-0.0044,-0.0051,-0.0060,-0.0068,-0.0078,-0.0088,-0.0098,-0.0107,-0.0116,-0.0124,-0.0132,-0.0138,-0.0144,-0.0149,-0.0154,-0.0160,-0.0166,-0.0172,-0.0178,-0.0185,-0.0191,-0.0196,-0.0201,
                                                                0.0501,0.0447,0.0383,0.0315 ,0.0248 ,0.0188 ,0.0137 ,0.0097 ,0.0066 ,0.0045 ,0.0029 ,0.0018 ,0.0011 ,0.0006 ,0.0002 ,-0.0001,-0.0003,-0.0005,-0.0007,-0.0009,-0.0011,-0.0013,-0.0016,-0.0019,-0.0023,-0.0027,-0.0032,-0.0038,-0.0045,-0.0052,-0.0060,-0.0069,-0.0078,-0.0088,-0.0097,-0.0106,-0.0115,-0.0122,-0.0129,-0.0135,-0.0141,-0.0146,-0.0152,-0.0158,-0.0165,-0.0171,-0.0178,-0.0185,-0.0191,-0.0197,
                                                                0.0518,0.0467,0.0406,0.0339 ,0.0271 ,0.0208 ,0.0153 ,0.0109 ,0.0076 ,0.0051 ,0.0034 ,0.0022 ,0.0014 ,0.0008 ,0.0004 ,0.0001 ,-0.0001,-0.0003,-0.0005,-0.0007,-0.0009,-0.0011,-0.0013,-0.0016,-0.0019,-0.0023,-0.0028,-0.0033,-0.0039,-0.0045,-0.0053,-0.0061,-0.0069,-0.0078,-0.0088,-0.0096,-0.0105,-0.0113,-0.0120,-0.0126,-0.0132,-0.0138,-0.0143,-0.0150,-0.0157,-0.0164,-0.0171,-0.0178,-0.0185,-0.0191,
                                                                0.0533,0.0485,0.0428,0.0363 ,0.0295 ,0.0229 ,0.0171 ,0.0124 ,0.0087 ,0.0059 ,0.0040 ,0.0026 ,0.0017 ,0.0010 ,0.0006 ,0.0002 ,-0.0000,-0.0002,-0.0004,-0.0005,-0.0007,-0.0009,-0.0011,-0.0013,-0.0016,-0.0020,-0.0024,-0.0028,-0.0033,-0.0039,-0.0046,-0.0053,-0.0061,-0.0070,-0.0078,-0.0087,-0.0095,-0.0103,-0.0110,-0.0117,-0.0123,-0.0129,-0.0135,-0.0141,-0.0148,-0.0155,-0.0163,-0.0171,-0.0178,-0.0185,
                                                                0.0547,0.0503,0.0449,0.0386 ,0.0318 ,0.0251 ,0.0190 ,0.0139 ,0.0098 ,0.0068 ,0.0046 ,0.0031 ,0.0020 ,0.0013 ,0.0007 ,0.0004 ,0.0001 ,-0.0001,-0.0002,-0.0004,-0.0005,-0.0007,-0.0009,-0.0011,-0.0014,-0.0017,-0.0020,-0.0024,-0.0029,-0.0034,-0.0040,-0.0047,-0.0054,-0.0062,-0.0070,-0.0078,-0.0086,-0.0094,-0.0101,-0.0108,-0.0114,-0.0120,-0.0126,-0.0132,-0.0139,-0.0147,-0.0155,-0.0163,-0.0171,-0.0178,
                                                                0.0559,0.0519,0.0469,0.0409 ,0.0342 ,0.0274 ,0.0211 ,0.0156 ,0.0112 ,0.0078 ,0.0053 ,0.0036 ,0.0024 ,0.0015 ,0.0010 ,0.0006 ,0.0003 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0006,-0.0007,-0.0009,-0.0011,-0.0014,-0.0017,-0.0020,-0.0024,-0.0029,-0.0034,-0.0040,-0.0047,-0.0054,-0.0062,-0.0070,-0.0077,-0.0085,-0.0092,-0.0099,-0.0105,-0.0111,-0.0117,-0.0123,-0.0130,-0.0138,-0.0146,-0.0154,-0.0163,-0.0171,
                                                                0.0571,0.0534,0.0488,0.0431 ,0.0366 ,0.0297 ,0.0232 ,0.0174 ,0.0126 ,0.0089 ,0.0061 ,0.0041 ,0.0028 ,0.0018 ,0.0012 ,0.0007 ,0.0004 ,0.0002 ,0.0000 ,-0.0001,-0.0003,-0.0004,-0.0006,-0.0007,-0.0009,-0.0012,-0.0014,-0.0017,-0.0021,-0.0025,-0.0030,-0.0035,-0.0041,-0.0047,-0.0054,-0.0062,-0.0069,-0.0076,-0.0083,-0.0090,-0.0096,-0.0102,-0.0108,-0.0114,-0.0121,-0.0128,-0.0137,-0.0145,-0.0154,-0.0163,
                                                                0.0581,0.0548,0.0505,0.0452 ,0.0389 ,0.0321 ,0.0254 ,0.0193 ,0.0141 ,0.0101 ,0.0070 ,0.0048 ,0.0032 ,0.0021 ,0.0014 ,0.0009 ,0.0006 ,0.0003 ,0.0001 ,-0.0000,-0.0002,-0.0003,-0.0004,-0.0006,-0.0007,-0.0009,-0.0012,-0.0014,-0.0018,-0.0021,-0.0025,-0.0030,-0.0035,-0.0041,-0.0048,-0.0054,-0.0061,-0.0068,-0.0075,-0.0081,-0.0087,-0.0093,-0.0099,-0.0105,-0.0111,-0.0119,-0.0127,-0.0136,-0.0145,-0.0154,
                                                                0.0590,0.0561,0.0522,0.0471 ,0.0412 ,0.0345 ,0.0277 ,0.0214 ,0.0158 ,0.0114 ,0.0080 ,0.0055 ,0.0037 ,0.0025 ,0.0017 ,0.0011 ,0.0007 ,0.0004 ,0.0002 ,0.0001 ,-0.0001,-0.0002,-0.0003,-0.0004,-0.0006,-0.0008,-0.0010,-0.0012,-0.0015,-0.0018,-0.0022,-0.0026,-0.0031,-0.0036,-0.0042,-0.0048,-0.0054,-0.0061,-0.0067,-0.0073,-0.0079,-0.0085,-0.0090,-0.0096,-0.0102,-0.0109,-0.0117,-0.0126,-0.0135,-0.0145,
                                                                0.0598,0.0572,0.0537,0.0490 ,0.0433 ,0.0369 ,0.0301 ,0.0235 ,0.0177 ,0.0128 ,0.0091 ,0.0063 ,0.0043 ,0.0029 ,0.0020 ,0.0013 ,0.0009 ,0.0006 ,0.0003 ,0.0002 ,0.0000 ,-0.0001,-0.0002,-0.0003,-0.0005,-0.0006,-0.0008,-0.0010,-0.0012,-0.0015,-0.0018,-0.0022,-0.0026,-0.0031,-0.0036,-0.0042,-0.0048,-0.0054,-0.0060,-0.0066,-0.0071,-0.0076,-0.0082,-0.0087,-0.0093,-0.0099,-0.0107,-0.0116,-0.0125,-0.0135,
                                                                0.0606,0.0582,0.0550,0.0508 ,0.0454 ,0.0392 ,0.0325 ,0.0257 ,0.0196 ,0.0144 ,0.0103 ,0.0072 ,0.0050 ,0.0034 ,0.0023 ,0.0015 ,0.0010 ,0.0007 ,0.0004 ,0.0003 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0003,-0.0005,-0.0006,-0.0008,-0.0010,-0.0013,-0.0015,-0.0019,-0.0022,-0.0027,-0.0031,-0.0036,-0.0042,-0.0047,-0.0053,-0.0058,-0.0064,-0.0069,-0.0074,-0.0079,-0.0084,-0.0090,-0.0097,-0.0106,-0.0115,-0.0125,
                                                                0.0612,0.0592,0.0563,0.0524 ,0.0474 ,0.0415 ,0.0348 ,0.0280 ,0.0217 ,0.0161 ,0.0116 ,0.0082 ,0.0057 ,0.0039 ,0.0026 ,0.0018 ,0.0012 ,0.0008 ,0.0005 ,0.0004 ,0.0002 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0004,-0.0005,-0.0007,-0.0008,-0.0010,-0.0013,-0.0016,-0.0019,-0.0023,-0.0027,-0.0031,-0.0036,-0.0041,-0.0046,-0.0052,-0.0057,-0.0061,-0.0066,-0.0070,-0.0075,-0.0081,-0.0088,-0.0095,-0.0104,-0.0114,
                                                                0.0618,0.0600,0.0574,0.0539 ,0.0493 ,0.0436 ,0.0372 ,0.0304 ,0.0238 ,0.0179 ,0.0131 ,0.0093 ,0.0065 ,0.0045 ,0.0031 ,0.0021 ,0.0014 ,0.0010 ,0.0007 ,0.0004 ,0.0003 ,0.0002 ,0.0000 ,-0.0001,-0.0002,-0.0003,-0.0004,-0.0005,-0.0007,-0.0009,-0.0011,-0.0013,-0.0016,-0.0019,-0.0023,-0.0027,-0.0031,-0.0036,-0.0041,-0.0045,-0.0050,-0.0055,-0.0059,-0.0063,-0.0067,-0.0072,-0.0078,-0.0085,-0.0094,-0.0104,
                                                                0.0623,0.0607,0.0584,0.0552 ,0.0510 ,0.0457 ,0.0395 ,0.0328 ,0.0261 ,0.0199 ,0.0147 ,0.0105 ,0.0074 ,0.0051 ,0.0035 ,0.0024 ,0.0016 ,0.0011 ,0.0008 ,0.0005 ,0.0004 ,0.0002 ,0.0001 ,0.0000 ,-0.0001,-0.0002,-0.0003,-0.0004,-0.0005,-0.0007,-0.0009,-0.0011,-0.0013,-0.0016,-0.0020,-0.0023,-0.0027,-0.0031,-0.0035,-0.0040,-0.0044,-0.0048,-0.0052,-0.0056,-0.0060,-0.0064,-0.0069,-0.0076,-0.0083,-0.0093,
                                                                0.0628,0.0614,0.0593,0.0565 ,0.0526 ,0.0477 ,0.0418 ,0.0352 ,0.0284 ,0.0220 ,0.0164 ,0.0118 ,0.0084 ,0.0058 ,0.0040 ,0.0028 ,0.0019 ,0.0013 ,0.0009 ,0.0006 ,0.0004 ,0.0003 ,0.0002 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0003,-0.0004,-0.0006,-0.0007,-0.0009,-0.0011,-0.0014,-0.0017,-0.0020,-0.0023,-0.0027,-0.0031,-0.0035,-0.0039,-0.0042,-0.0046,-0.0049,-0.0053,-0.0056,-0.0061,-0.0066,-0.0073,-0.0082,
                                                                0.0632,0.0620,0.0602,0.0576 ,0.0541 ,0.0496 ,0.0439 ,0.0375 ,0.0307 ,0.0241 ,0.0182 ,0.0133 ,0.0095 ,0.0066 ,0.0046 ,0.0032 ,0.0022 ,0.0015 ,0.0010 ,0.0007 ,0.0005 ,0.0004 ,0.0002 ,0.0001 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0003,-0.0004,-0.0006,-0.0007,-0.0009,-0.0011,-0.0014,-0.0017,-0.0020,-0.0023,-0.0026,-0.0030,-0.0034,-0.0037,-0.0040,-0.0043,-0.0046,-0.0049,-0.0053,-0.0058,-0.0064,-0.0072};
        
        return arrayCurvatureMin[y * arrayWidth + x];
    }
    case 7:
    {
        static const std::vector<double> arrayCurvatureMin   = {0.0030,0.0020,0.0012,0.0006,0.0000,-0.0005,-0.0011,-0.0018,-0.0026,-0.0037,-0.0049,-0.0062,-0.0076,-0.0090,-0.0103,-0.0115,-0.0126,-0.0135,-0.0142,-0.0148,-0.0152,-0.0156,-0.0158,-0.0160,-0.0161,-0.0162,-0.0163,-0.0163,-0.0163,-0.0163,-0.0163,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0162,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0160,-0.0160,-0.0160,-0.0160,-0.0160,
                                                                0.0043,0.0031,0.0021,0.0013,0.0006,0.0000 ,-0.0006,-0.0013,-0.0022,-0.0032,-0.0043,-0.0056,-0.0070,-0.0084,-0.0098,-0.0111,-0.0123,-0.0133,-0.0141,-0.0147,-0.0152,-0.0156,-0.0159,-0.0161,-0.0162,-0.0163,-0.0164,-0.0164,-0.0164,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0162,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0160,-0.0160,-0.0160,-0.0160,
                                                                0.0056,0.0041,0.0029,0.0019,0.0011,0.0004 ,-0.0004,-0.0011,-0.0019,-0.0028,-0.0039,-0.0051,-0.0064,-0.0077,-0.0091,-0.0105,-0.0117,-0.0128,-0.0138,-0.0145,-0.0151,-0.0156,-0.0159,-0.0161,-0.0163,-0.0164,-0.0164,-0.0165,-0.0165,-0.0165,-0.0164,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0160,-0.0160,-0.0160,-0.0160,
                                                                0.0068,0.0051,0.0036,0.0024,0.0014,0.0006 ,-0.0002,-0.0010,-0.0018,-0.0026,-0.0035,-0.0046,-0.0058,-0.0070,-0.0084,-0.0097,-0.0110,-0.0122,-0.0132,-0.0141,-0.0148,-0.0153,-0.0157,-0.0160,-0.0162,-0.0164,-0.0165,-0.0165,-0.0165,-0.0165,-0.0165,-0.0165,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0160,-0.0160,-0.0160,
                                                                0.0079,0.0059,0.0043,0.0029,0.0017,0.0007 ,-0.0002,-0.0010,-0.0017,-0.0025,-0.0033,-0.0042,-0.0052,-0.0063,-0.0075,-0.0088,-0.0101,-0.0113,-0.0124,-0.0134,-0.0142,-0.0149,-0.0154,-0.0158,-0.0161,-0.0163,-0.0164,-0.0165,-0.0165,-0.0165,-0.0165,-0.0165,-0.0165,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0162,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0160,-0.0160,
                                                                0.0090,0.0068,0.0048,0.0032,0.0019,0.0007 ,-0.0002,-0.0011,-0.0018,-0.0025,-0.0032,-0.0039,-0.0048,-0.0057,-0.0067,-0.0078,-0.0090,-0.0102,-0.0114,-0.0125,-0.0135,-0.0143,-0.0149,-0.0155,-0.0158,-0.0161,-0.0163,-0.0164,-0.0165,-0.0165,-0.0165,-0.0165,-0.0165,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0162,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0160,
                                                                0.0100,0.0075,0.0054,0.0036,0.0020,0.0008 ,-0.0003,-0.0011,-0.0019,-0.0025,-0.0031,-0.0037,-0.0044,-0.0051,-0.0059,-0.0069,-0.0080,-0.0091,-0.0103,-0.0115,-0.0125,-0.0135,-0.0143,-0.0149,-0.0154,-0.0158,-0.0161,-0.0163,-0.0164,-0.0165,-0.0165,-0.0165,-0.0165,-0.0165,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0162,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,
                                                                0.0109,0.0083,0.0059,0.0039,0.0022,0.0008 ,-0.0003,-0.0012,-0.0019,-0.0025,-0.0030,-0.0035,-0.0040,-0.0046,-0.0053,-0.0061,-0.0070,-0.0080,-0.0091,-0.0103,-0.0114,-0.0125,-0.0134,-0.0142,-0.0149,-0.0154,-0.0158,-0.0160,-0.0162,-0.0164,-0.0164,-0.0165,-0.0165,-0.0165,-0.0165,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0162,-0.0161,-0.0161,-0.0161,-0.0161,-0.0161,
                                                                0.0118,0.0090,0.0064,0.0042,0.0024,0.0009 ,-0.0003,-0.0012,-0.0019,-0.0025,-0.0030,-0.0034,-0.0038,-0.0042,-0.0047,-0.0053,-0.0060,-0.0069,-0.0079,-0.0090,-0.0102,-0.0113,-0.0124,-0.0133,-0.0141,-0.0148,-0.0153,-0.0157,-0.0160,-0.0162,-0.0163,-0.0164,-0.0164,-0.0165,-0.0165,-0.0164,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0162,-0.0161,-0.0161,-0.0161,-0.0161,
                                                                0.0128,0.0097,0.0070,0.0046,0.0027,0.0011 ,-0.0002,-0.0012,-0.0019,-0.0025,-0.0029,-0.0032,-0.0035,-0.0038,-0.0042,-0.0047,-0.0052,-0.0059,-0.0068,-0.0078,-0.0089,-0.0100,-0.0112,-0.0122,-0.0132,-0.0140,-0.0147,-0.0152,-0.0156,-0.0159,-0.0162,-0.0163,-0.0164,-0.0164,-0.0164,-0.0165,-0.0164,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0162,-0.0161,-0.0161,-0.0161,
                                                                0.0137,0.0105,0.0076,0.0051,0.0030,0.0013 ,-0.0000,-0.0010,-0.0018,-0.0024,-0.0028,-0.0031,-0.0033,-0.0035,-0.0038,-0.0041,-0.0045,-0.0051,-0.0058,-0.0066,-0.0076,-0.0087,-0.0099,-0.0110,-0.0121,-0.0131,-0.0139,-0.0146,-0.0152,-0.0156,-0.0159,-0.0161,-0.0163,-0.0164,-0.0164,-0.0164,-0.0164,-0.0164,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0162,-0.0161,-0.0161,
                                                                0.0146,0.0113,0.0083,0.0056,0.0034,0.0017 ,0.0002 ,-0.0008,-0.0016,-0.0022,-0.0026,-0.0029,-0.0031,-0.0033,-0.0034,-0.0037,-0.0039,-0.0044,-0.0049,-0.0056,-0.0064,-0.0074,-0.0086,-0.0097,-0.0109,-0.0120,-0.0130,-0.0138,-0.0145,-0.0151,-0.0155,-0.0159,-0.0161,-0.0162,-0.0163,-0.0164,-0.0164,-0.0164,-0.0164,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,-0.0162,
                                                                0.0156,0.0122,0.0090,0.0063,0.0040,0.0021 ,0.0006 ,-0.0005,-0.0014,-0.0020,-0.0024,-0.0027,-0.0029,-0.0030,-0.0031,-0.0033,-0.0035,-0.0037,-0.0041,-0.0047,-0.0054,-0.0062,-0.0073,-0.0084,-0.0096,-0.0107,-0.0119,-0.0129,-0.0138,-0.0145,-0.0151,-0.0155,-0.0158,-0.0161,-0.0162,-0.0163,-0.0164,-0.0164,-0.0165,-0.0165,-0.0164,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0162,-0.0162,-0.0162,
                                                                0.0167,0.0131,0.0099,0.0070,0.0046,0.0026 ,0.0010 ,-0.0002,-0.0011,-0.0017,-0.0022,-0.0025,-0.0026,-0.0028,-0.0029,-0.0029,-0.0031,-0.0033,-0.0035,-0.0039,-0.0045,-0.0052,-0.0061,-0.0071,-0.0082,-0.0094,-0.0106,-0.0117,-0.0128,-0.0137,-0.0144,-0.0150,-0.0155,-0.0158,-0.0160,-0.0162,-0.0163,-0.0164,-0.0165,-0.0165,-0.0165,-0.0165,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,-0.0163,-0.0162,
                                                                0.0178,0.0142,0.0108,0.0079,0.0053,0.0033 ,0.0016 ,0.0003 ,-0.0007,-0.0014,-0.0019,-0.0022,-0.0024,-0.0025,-0.0026,-0.0027,-0.0027,-0.0028,-0.0030,-0.0033,-0.0037,-0.0043,-0.0050,-0.0059,-0.0069,-0.0081,-0.0093,-0.0105,-0.0116,-0.0127,-0.0136,-0.0144,-0.0150,-0.0154,-0.0158,-0.0160,-0.0162,-0.0163,-0.0164,-0.0165,-0.0165,-0.0165,-0.0165,-0.0165,-0.0164,-0.0164,-0.0164,-0.0163,-0.0163,-0.0163,
                                                                0.0191,0.0153,0.0119,0.0088,0.0062,0.0040 ,0.0022 ,0.0008 ,-0.0002,-0.0010,-0.0015,-0.0019,-0.0021,-0.0023,-0.0023,-0.0024,-0.0024,-0.0025,-0.0026,-0.0028,-0.0031,-0.0035,-0.0041,-0.0048,-0.0057,-0.0068,-0.0079,-0.0091,-0.0104,-0.0115,-0.0126,-0.0135,-0.0143,-0.0149,-0.0154,-0.0158,-0.0160,-0.0162,-0.0163,-0.0164,-0.0165,-0.0165,-0.0165,-0.0165,-0.0165,-0.0164,-0.0164,-0.0164,-0.0164,-0.0163,
                                                                0.0204,0.0166,0.0131,0.0099,0.0071,0.0048 ,0.0029 ,0.0014 ,0.0003 ,-0.0006,-0.0012,-0.0016,-0.0019,-0.0020,-0.0021,-0.0021,-0.0022,-0.0022,-0.0023,-0.0024,-0.0026,-0.0029,-0.0033,-0.0039,-0.0047,-0.0056,-0.0066,-0.0078,-0.0090,-0.0103,-0.0115,-0.0125,-0.0135,-0.0143,-0.0149,-0.0154,-0.0158,-0.0160,-0.0162,-0.0164,-0.0164,-0.0165,-0.0165,-0.0165,-0.0165,-0.0165,-0.0165,-0.0164,-0.0164,-0.0164,
                                                                0.0218,0.0179,0.0143,0.0110,0.0082,0.0057 ,0.0037 ,0.0021 ,0.0009 ,-0.0001,-0.0008,-0.0012,-0.0016,-0.0017,-0.0019,-0.0019,-0.0019,-0.0020,-0.0020,-0.0021,-0.0022,-0.0024,-0.0027,-0.0032,-0.0038,-0.0045,-0.0054,-0.0065,-0.0077,-0.0089,-0.0102,-0.0114,-0.0125,-0.0134,-0.0142,-0.0149,-0.0154,-0.0157,-0.0160,-0.0162,-0.0164,-0.0165,-0.0165,-0.0165,-0.0166,-0.0165,-0.0165,-0.0165,-0.0165,-0.0164,
                                                                0.0233,0.0194,0.0157,0.0123,0.0093,0.0067 ,0.0046 ,0.0028 ,0.0015 ,0.0004 ,-0.0003,-0.0009,-0.0012,-0.0015,-0.0016,-0.0017,-0.0017,-0.0017,-0.0018,-0.0018,-0.0019,-0.0020,-0.0022,-0.0026,-0.0030,-0.0036,-0.0044,-0.0053,-0.0064,-0.0076,-0.0088,-0.0101,-0.0113,-0.0124,-0.0134,-0.0142,-0.0148,-0.0153,-0.0157,-0.0160,-0.0162,-0.0164,-0.0165,-0.0165,-0.0166,-0.0166,-0.0166,-0.0166,-0.0165,-0.0165,
                                                                0.0248,0.0209,0.0171,0.0136,0.0105,0.0078 ,0.0055 ,0.0036 ,0.0022 ,0.0010 ,0.0002 ,-0.0005,-0.0009,-0.0012,-0.0014,-0.0015,-0.0015,-0.0015,-0.0016,-0.0016,-0.0016,-0.0017,-0.0019,-0.0021,-0.0024,-0.0029,-0.0035,-0.0043,-0.0052,-0.0063,-0.0075,-0.0087,-0.0100,-0.0112,-0.0123,-0.0133,-0.0141,-0.0148,-0.0153,-0.0157,-0.0160,-0.0163,-0.0164,-0.0165,-0.0166,-0.0166,-0.0166,-0.0166,-0.0166,-0.0166,
                                                                0.0265,0.0225,0.0186,0.0150,0.0117,0.0089 ,0.0065 ,0.0045 ,0.0029 ,0.0016 ,0.0007 ,-0.0000,-0.0005,-0.0009,-0.0011,-0.0012,-0.0013,-0.0013,-0.0014,-0.0014,-0.0014,-0.0015,-0.0016,-0.0017,-0.0020,-0.0023,-0.0028,-0.0034,-0.0042,-0.0051,-0.0062,-0.0074,-0.0087,-0.0099,-0.0112,-0.0123,-0.0133,-0.0141,-0.0148,-0.0153,-0.0157,-0.0160,-0.0163,-0.0164,-0.0165,-0.0166,-0.0166,-0.0166,-0.0166,-0.0166,
                                                                0.0282,0.0242,0.0202,0.0165,0.0131,0.0101 ,0.0075 ,0.0054 ,0.0036 ,0.0023 ,0.0012 ,0.0004 ,-0.0001,-0.0006,-0.0008,-0.0010,-0.0011,-0.0012,-0.0012,-0.0012,-0.0012,-0.0012,-0.0013,-0.0014,-0.0016,-0.0018,-0.0022,-0.0027,-0.0033,-0.0041,-0.0050,-0.0061,-0.0073,-0.0086,-0.0099,-0.0111,-0.0122,-0.0132,-0.0141,-0.0148,-0.0153,-0.0157,-0.0161,-0.0163,-0.0164,-0.0166,-0.0166,-0.0167,-0.0167,-0.0167,
                                                                0.0299,0.0259,0.0219,0.0180,0.0145,0.0113 ,0.0086 ,0.0063 ,0.0044 ,0.0029 ,0.0018 ,0.0009 ,0.0003 ,-0.0002,-0.0005,-0.0008,-0.0009,-0.0010,-0.0010,-0.0010,-0.0011,-0.0011,-0.0011,-0.0012,-0.0013,-0.0015,-0.0017,-0.0021,-0.0026,-0.0032,-0.0040,-0.0049,-0.0060,-0.0072,-0.0085,-0.0098,-0.0111,-0.0122,-0.0132,-0.0141,-0.0148,-0.0153,-0.0157,-0.0161,-0.0163,-0.0165,-0.0166,-0.0166,-0.0167,-0.0167,
                                                                0.0317,0.0277,0.0236,0.0196,0.0159,0.0126 ,0.0097 ,0.0073 ,0.0053 ,0.0037 ,0.0024 ,0.0014 ,0.0007 ,0.0001 ,-0.0003,-0.0005,-0.0007,-0.0008,-0.0009,-0.0009,-0.0009,-0.0009,-0.0009,-0.0010,-0.0011,-0.0012,-0.0014,-0.0016,-0.0020,-0.0025,-0.0031,-0.0039,-0.0048,-0.0059,-0.0072,-0.0084,-0.0097,-0.0110,-0.0122,-0.0132,-0.0140,-0.0147,-0.0153,-0.0157,-0.0161,-0.0163,-0.0165,-0.0166,-0.0167,-0.0167,
                                                                0.0336,0.0295,0.0253,0.0212,0.0174,0.0140 ,0.0109 ,0.0083 ,0.0062 ,0.0044 ,0.0030 ,0.0019 ,0.0011 ,0.0005 ,0.0000 ,-0.0003,-0.0005,-0.0006,-0.0007,-0.0008,-0.0008,-0.0008,-0.0008,-0.0008,-0.0009,-0.0010,-0.0011,-0.0013,-0.0016,-0.0019,-0.0024,-0.0030,-0.0038,-0.0048,-0.0059,-0.0071,-0.0084,-0.0097,-0.0110,-0.0121,-0.0131,-0.0140,-0.0147,-0.0153,-0.0157,-0.0161,-0.0163,-0.0165,-0.0166,-0.0167,
                                                                0.0354,0.0313,0.0271,0.0229,0.0190,0.0153 ,0.0121 ,0.0094 ,0.0071 ,0.0052 ,0.0037 ,0.0025 ,0.0016 ,0.0009 ,0.0003 ,-0.0000,-0.0003,-0.0004,-0.0006,-0.0006,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0008,-0.0009,-0.0010,-0.0012,-0.0015,-0.0019,-0.0023,-0.0030,-0.0038,-0.0047,-0.0058,-0.0070,-0.0083,-0.0096,-0.0109,-0.0121,-0.0131,-0.0140,-0.0147,-0.0153,-0.0157,-0.0161,-0.0163,-0.0165,-0.0166,
                                                                0.0373,0.0332,0.0289,0.0246,0.0205,0.0168 ,0.0134 ,0.0105 ,0.0080 ,0.0060 ,0.0043 ,0.0030 ,0.0020 ,0.0012 ,0.0007 ,0.0002 ,-0.0001,-0.0003,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0009,-0.0011,-0.0014,-0.0018,-0.0023,-0.0029,-0.0037,-0.0047,-0.0057,-0.0070,-0.0083,-0.0096,-0.0109,-0.0120,-0.0131,-0.0140,-0.0147,-0.0153,-0.0157,-0.0161,-0.0163,-0.0165,
                                                                0.0392,0.0351,0.0307,0.0264,0.0222,0.0182 ,0.0147 ,0.0116 ,0.0090 ,0.0068 ,0.0050 ,0.0036 ,0.0025 ,0.0016 ,0.0010 ,0.0005 ,0.0002 ,-0.0001,-0.0002,-0.0004,-0.0004,-0.0005,-0.0005,-0.0005,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0009,-0.0011,-0.0014,-0.0017,-0.0022,-0.0029,-0.0037,-0.0046,-0.0057,-0.0069,-0.0082,-0.0095,-0.0108,-0.0120,-0.0131,-0.0140,-0.0147,-0.0153,-0.0157,-0.0161,-0.0164,
                                                                0.0411,0.0370,0.0326,0.0281,0.0238,0.0197 ,0.0160 ,0.0128 ,0.0100 ,0.0077 ,0.0058 ,0.0042 ,0.0030 ,0.0021 ,0.0013 ,0.0008 ,0.0004 ,0.0001 ,-0.0001,-0.0002,-0.0003,-0.0004,-0.0004,-0.0004,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0007,-0.0008,-0.0010,-0.0013,-0.0017,-0.0022,-0.0028,-0.0036,-0.0045,-0.0056,-0.0069,-0.0082,-0.0095,-0.0108,-0.0120,-0.0130,-0.0139,-0.0147,-0.0153,-0.0157,-0.0161,
                                                                0.0430,0.0389,0.0345,0.0299,0.0255,0.0213 ,0.0174 ,0.0140 ,0.0111 ,0.0086 ,0.0065 ,0.0049 ,0.0035 ,0.0025 ,0.0017 ,0.0011 ,0.0006 ,0.0003 ,0.0001 ,-0.0001,-0.0002,-0.0003,-0.0003,-0.0003,-0.0004,-0.0004,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0008,-0.0010,-0.0013,-0.0016,-0.0021,-0.0028,-0.0036,-0.0045,-0.0056,-0.0068,-0.0081,-0.0094,-0.0107,-0.0119,-0.0130,-0.0139,-0.0147,-0.0153,-0.0157,
                                                                0.0449,0.0408,0.0364,0.0317,0.0272,0.0228 ,0.0188 ,0.0152 ,0.0121 ,0.0095 ,0.0073 ,0.0055 ,0.0041 ,0.0030 ,0.0021 ,0.0014 ,0.0009 ,0.0005 ,0.0003 ,0.0001 ,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0006,-0.0007,-0.0009,-0.0012,-0.0016,-0.0021,-0.0027,-0.0035,-0.0045,-0.0055,-0.0068,-0.0081,-0.0094,-0.0107,-0.0119,-0.0130,-0.0139,-0.0146,-0.0153,
                                                                0.0468,0.0427,0.0382,0.0336,0.0289,0.0244 ,0.0203 ,0.0165 ,0.0132 ,0.0104 ,0.0081 ,0.0062 ,0.0046 ,0.0034 ,0.0025 ,0.0017 ,0.0012 ,0.0007 ,0.0004 ,0.0002 ,0.0000 ,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0007,-0.0009,-0.0012,-0.0016,-0.0021,-0.0027,-0.0035,-0.0044,-0.0055,-0.0067,-0.0080,-0.0094,-0.0107,-0.0119,-0.0129,-0.0139,-0.0146,
                                                                0.0487,0.0446,0.0401,0.0354,0.0307,0.0260 ,0.0217 ,0.0178 ,0.0144 ,0.0114 ,0.0089 ,0.0069 ,0.0052 ,0.0039 ,0.0029 ,0.0020 ,0.0014 ,0.0010 ,0.0006 ,0.0004 ,0.0002 ,0.0000 ,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0005,-0.0007,-0.0009,-0.0012,-0.0015,-0.0020,-0.0027,-0.0034,-0.0044,-0.0055,-0.0067,-0.0080,-0.0093,-0.0106,-0.0118,-0.0129,-0.0138,
                                                                0.0506,0.0465,0.0420,0.0373,0.0324,0.0277 ,0.0232 ,0.0192 ,0.0156 ,0.0124 ,0.0098 ,0.0076 ,0.0058 ,0.0044 ,0.0033 ,0.0024 ,0.0017 ,0.0012 ,0.0008 ,0.0005 ,0.0003 ,0.0001 ,0.0000 ,-0.0000,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0002,-0.0002,-0.0003,-0.0004,-0.0005,-0.0006,-0.0009,-0.0011,-0.0015,-0.0020,-0.0026,-0.0034,-0.0043,-0.0054,-0.0066,-0.0079,-0.0093,-0.0106,-0.0118,-0.0129,
                                                                0.0524,0.0484,0.0439,0.0391,0.0342,0.0293 ,0.0248 ,0.0205 ,0.0168 ,0.0135 ,0.0107 ,0.0084 ,0.0065 ,0.0049 ,0.0037 ,0.0027 ,0.0020 ,0.0014 ,0.0010 ,0.0007 ,0.0004 ,0.0002 ,0.0001 ,0.0000 ,-0.0000,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0005,-0.0006,-0.0008,-0.0011,-0.0015,-0.0020,-0.0026,-0.0034,-0.0043,-0.0054,-0.0066,-0.0079,-0.0092,-0.0105,-0.0118,
                                                                0.0542,0.0503,0.0458,0.0410,0.0360,0.0310 ,0.0263 ,0.0219 ,0.0180 ,0.0146 ,0.0116 ,0.0091 ,0.0071 ,0.0055 ,0.0041 ,0.0031 ,0.0023 ,0.0017 ,0.0012 ,0.0008 ,0.0006 ,0.0004 ,0.0002 ,0.0001 ,0.0000 ,-0.0000,-0.0000,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0003,-0.0004,-0.0006,-0.0008,-0.0011,-0.0015,-0.0019,-0.0026,-0.0033,-0.0043,-0.0053,-0.0065,-0.0078,-0.0092,-0.0105,
                                                                0.0560,0.0521,0.0477,0.0428,0.0378,0.0327 ,0.0279 ,0.0234 ,0.0193 ,0.0157 ,0.0126 ,0.0099 ,0.0078 ,0.0060 ,0.0046 ,0.0035 ,0.0026 ,0.0019 ,0.0014 ,0.0010 ,0.0007 ,0.0005 ,0.0003 ,0.0002 ,0.0001 ,0.0000 ,0.0000 ,-0.0000,-0.0000,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0003,-0.0004,-0.0006,-0.0008,-0.0011,-0.0014,-0.0019,-0.0025,-0.0033,-0.0042,-0.0053,-0.0065,-0.0078,-0.0091,
                                                                0.0577,0.0539,0.0495,0.0447,0.0396,0.0345 ,0.0295 ,0.0248 ,0.0206 ,0.0168 ,0.0135 ,0.0108 ,0.0085 ,0.0066 ,0.0051 ,0.0039 ,0.0029 ,0.0022 ,0.0016 ,0.0012 ,0.0008 ,0.0006 ,0.0004 ,0.0003 ,0.0002 ,0.0001 ,0.0000 ,0.0000 ,-0.0000,-0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0001,-0.0001,-0.0002,-0.0003,-0.0004,-0.0006,-0.0008,-0.0010,-0.0014,-0.0019,-0.0025,-0.0033,-0.0042,-0.0053,-0.0065,-0.0078,
                                                                0.0594,0.0557,0.0513,0.0465,0.0414,0.0362 ,0.0311 ,0.0263 ,0.0219 ,0.0180 ,0.0145 ,0.0116 ,0.0092 ,0.0072 ,0.0056 ,0.0043 ,0.0032 ,0.0024 ,0.0018 ,0.0013 ,0.0010 ,0.0007 ,0.0005 ,0.0003 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,0.0000 ,-0.0000,-0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0001,-0.0002,-0.0003,-0.0004,-0.0005,-0.0007,-0.0010,-0.0014,-0.0019,-0.0025,-0.0032,-0.0041,-0.0052,-0.0064,
                                                                0.0611,0.0575,0.0532,0.0483,0.0432,0.0379 ,0.0327 ,0.0278 ,0.0232 ,0.0191 ,0.0156 ,0.0125 ,0.0099 ,0.0078 ,0.0061 ,0.0047 ,0.0036 ,0.0027 ,0.0020 ,0.0015 ,0.0011 ,0.0008 ,0.0006 ,0.0004 ,0.0003 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,-0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0002,-0.0003,-0.0004,-0.0005,-0.0007,-0.0010,-0.0014,-0.0019,-0.0025,-0.0032,-0.0041,-0.0052,
                                                                0.0628,0.0592,0.0550,0.0502,0.0450,0.0397 ,0.0344 ,0.0293 ,0.0246 ,0.0204 ,0.0166 ,0.0134 ,0.0107 ,0.0085 ,0.0066 ,0.0051 ,0.0040 ,0.0030 ,0.0023 ,0.0017 ,0.0013 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,-0.0000,-0.0000,-0.0001,-0.0001,-0.0002,-0.0002,-0.0004,-0.0005,-0.0007,-0.0010,-0.0014,-0.0018,-0.0024,-0.0032,-0.0041,
                                                                0.0644,0.0609,0.0567,0.0520,0.0468,0.0414 ,0.0361 ,0.0309 ,0.0260 ,0.0216 ,0.0177 ,0.0144 ,0.0115 ,0.0091 ,0.0072 ,0.0056 ,0.0043 ,0.0033 ,0.0025 ,0.0019 ,0.0014 ,0.0011 ,0.0008 ,0.0006 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,0.0000 ,-0.0000,-0.0000,-0.0001,-0.0001,-0.0002,-0.0003,-0.0005,-0.0007,-0.0010,-0.0013,-0.0018,-0.0024,-0.0031,
                                                                0.0660,0.0626,0.0585,0.0537,0.0486,0.0432 ,0.0377 ,0.0325 ,0.0275 ,0.0229 ,0.0189 ,0.0153 ,0.0123 ,0.0098 ,0.0078 ,0.0061 ,0.0047 ,0.0036 ,0.0028 ,0.0021 ,0.0016 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,0.0000 ,0.0000 ,-0.0000,-0.0000,-0.0001,-0.0001,-0.0002,-0.0003,-0.0005,-0.0007,-0.0010,-0.0013,-0.0018,-0.0024,
                                                                0.0675,0.0642,0.0602,0.0555,0.0504,0.0449 ,0.0394 ,0.0340 ,0.0289 ,0.0242 ,0.0200 ,0.0163 ,0.0132 ,0.0105 ,0.0084 ,0.0066 ,0.0051 ,0.0040 ,0.0031 ,0.0023 ,0.0018 ,0.0014 ,0.0010 ,0.0008 ,0.0006 ,0.0004 ,0.0003 ,0.0003 ,0.0002 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,0.0000 ,0.0000 ,-0.0000,-0.0001,-0.0001,-0.0002,-0.0003,-0.0005,-0.0007,-0.0010,-0.0013,-0.0018,
                                                                0.0690,0.0658,0.0619,0.0573,0.0521,0.0467 ,0.0411 ,0.0357 ,0.0304 ,0.0256 ,0.0212 ,0.0174 ,0.0141 ,0.0113 ,0.0090 ,0.0071 ,0.0055 ,0.0043 ,0.0033 ,0.0026 ,0.0020 ,0.0015 ,0.0011 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,0.0000 ,0.0000 ,-0.0000,-0.0001,-0.0001,-0.0002,-0.0003,-0.0005,-0.0007,-0.0009,-0.0013,
                                                                0.0705,0.0674,0.0635,0.0590,0.0539,0.0485 ,0.0429 ,0.0373 ,0.0319 ,0.0270 ,0.0224 ,0.0184 ,0.0150 ,0.0121 ,0.0096 ,0.0076 ,0.0060 ,0.0047 ,0.0036 ,0.0028 ,0.0022 ,0.0017 ,0.0013 ,0.0010 ,0.0007 ,0.0006 ,0.0004 ,0.0003 ,0.0003 ,0.0002 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,0.0000 ,-0.0000,-0.0001,-0.0001,-0.0002,-0.0003,-0.0005,-0.0007,-0.0009,
                                                                0.0719,0.0689,0.0651,0.0607,0.0556,0.0502 ,0.0446 ,0.0389 ,0.0335 ,0.0284 ,0.0237 ,0.0195 ,0.0159 ,0.0129 ,0.0103 ,0.0082 ,0.0064 ,0.0051 ,0.0039 ,0.0031 ,0.0024 ,0.0018 ,0.0014 ,0.0011 ,0.0008 ,0.0006 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,0.0000 ,-0.0000,-0.0001,-0.0001,-0.0002,-0.0003,-0.0005,-0.0007,
                                                                0.0733,0.0704,0.0667,0.0623,0.0574,0.0520 ,0.0463 ,0.0406 ,0.0350 ,0.0298 ,0.0250 ,0.0207 ,0.0169 ,0.0137 ,0.0110 ,0.0088 ,0.0069 ,0.0054 ,0.0043 ,0.0033 ,0.0026 ,0.0020 ,0.0015 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,-0.0000,-0.0000,-0.0001,-0.0002,-0.0003,-0.0004,
                                                                0.0746,0.0718,0.0683,0.0640,0.0591,0.0537 ,0.0480 ,0.0423 ,0.0366 ,0.0312 ,0.0263 ,0.0218 ,0.0179 ,0.0146 ,0.0117 ,0.0094 ,0.0074 ,0.0059 ,0.0046 ,0.0036 ,0.0028 ,0.0022 ,0.0017 ,0.0013 ,0.0010 ,0.0008 ,0.0006 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,0.0000 ,-0.0000,-0.0001,-0.0002,-0.0003,
                                                                0.0759,0.0732,0.0698,0.0656,0.0608,0.0554 ,0.0497 ,0.0439 ,0.0382 ,0.0327 ,0.0276 ,0.0230 ,0.0190 ,0.0154 ,0.0125 ,0.0100 ,0.0079 ,0.0063 ,0.0049 ,0.0039 ,0.0030 ,0.0023 ,0.0018 ,0.0014 ,0.0011 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,0.0000 ,-0.0000,-0.0001,-0.0002};
        
        return arrayCurvatureMin[y * arrayWidth + x];
    }
    case 8:
    {
        static const std::vector<double> arrayCurvatureMin   = {0.0155,0.0110,0.0074,0.0047,0.0026,0.0012,0.0002,-0.0005,-0.0009,-0.0012,-0.0014,-0.0017,-0.0020,-0.0024,-0.0029,-0.0034,-0.0041,-0.0048,-0.0056,-0.0065,-0.0073,-0.0082,-0.0090,-0.0098,-0.0105,-0.0111,-0.0117,-0.0121,-0.0125,-0.0128,-0.0130,-0.0131,-0.0132,-0.0132,-0.0132,-0.0131,-0.0128,-0.0125,-0.0121,-0.0115,-0.0108,-0.0101,-0.0095,-0.0089,-0.0084,-0.0080,-0.0077,-0.0075,-0.0074,-0.0073,
                                                                0.0165,0.0118,0.0081,0.0052,0.0031,0.0015,0.0005,-0.0002,-0.0006,-0.0009,-0.0011,-0.0013,-0.0016,-0.0019,-0.0023,-0.0028,-0.0033,-0.0040,-0.0047,-0.0055,-0.0063,-0.0072,-0.0081,-0.0089,-0.0097,-0.0104,-0.0111,-0.0116,-0.0121,-0.0124,-0.0127,-0.0129,-0.0130,-0.0131,-0.0131,-0.0130,-0.0128,-0.0126,-0.0121,-0.0116,-0.0110,-0.0103,-0.0096,-0.0090,-0.0085,-0.0081,-0.0078,-0.0076,-0.0074,-0.0073,
                                                                0.0175,0.0126,0.0087,0.0057,0.0035,0.0019,0.0008,0.0001 ,-0.0003,-0.0006,-0.0008,-0.0010,-0.0013,-0.0015,-0.0018,-0.0022,-0.0027,-0.0033,-0.0039,-0.0046,-0.0054,-0.0062,-0.0071,-0.0080,-0.0088,-0.0096,-0.0103,-0.0110,-0.0115,-0.0120,-0.0123,-0.0126,-0.0128,-0.0130,-0.0130,-0.0130,-0.0128,-0.0126,-0.0122,-0.0117,-0.0111,-0.0104,-0.0097,-0.0091,-0.0086,-0.0081,-0.0078,-0.0076,-0.0074,-0.0073,
                                                                0.0186,0.0135,0.0095,0.0063,0.0040,0.0023,0.0011,0.0004 ,-0.0001,-0.0004,-0.0006,-0.0008,-0.0010,-0.0012,-0.0015,-0.0018,-0.0022,-0.0027,-0.0032,-0.0038,-0.0045,-0.0053,-0.0061,-0.0070,-0.0079,-0.0087,-0.0095,-0.0103,-0.0109,-0.0114,-0.0119,-0.0123,-0.0126,-0.0127,-0.0128,-0.0129,-0.0128,-0.0126,-0.0122,-0.0118,-0.0112,-0.0105,-0.0098,-0.0092,-0.0087,-0.0082,-0.0079,-0.0076,-0.0074,-0.0073,
                                                                0.0197,0.0144,0.0102,0.0069,0.0045,0.0027,0.0015,0.0007 ,0.0001 ,-0.0002,-0.0004,-0.0006,-0.0008,-0.0010,-0.0012,-0.0015,-0.0018,-0.0022,-0.0026,-0.0031,-0.0037,-0.0044,-0.0052,-0.0060,-0.0069,-0.0078,-0.0086,-0.0094,-0.0101,-0.0108,-0.0114,-0.0118,-0.0122,-0.0125,-0.0126,-0.0127,-0.0127,-0.0125,-0.0122,-0.0118,-0.0113,-0.0106,-0.0099,-0.0093,-0.0087,-0.0083,-0.0079,-0.0077,-0.0075,-0.0073,
                                                                0.0210,0.0155,0.0110,0.0076,0.0050,0.0031,0.0018,0.0009 ,0.0004 ,0.0000 ,-0.0002,-0.0004,-0.0006,-0.0008,-0.0010,-0.0012,-0.0014,-0.0018,-0.0021,-0.0026,-0.0031,-0.0037,-0.0043,-0.0051,-0.0059,-0.0068,-0.0076,-0.0085,-0.0093,-0.0100,-0.0107,-0.0113,-0.0117,-0.0121,-0.0123,-0.0125,-0.0125,-0.0124,-0.0122,-0.0118,-0.0113,-0.0107,-0.0100,-0.0094,-0.0088,-0.0083,-0.0080,-0.0077,-0.0075,-0.0074,
                                                                0.0223,0.0165,0.0119,0.0083,0.0056,0.0036,0.0022,0.0012 ,0.0006 ,0.0002 ,-0.0000,-0.0002,-0.0004,-0.0006,-0.0008,-0.0010,-0.0012,-0.0014,-0.0017,-0.0021,-0.0025,-0.0030,-0.0036,-0.0043,-0.0050,-0.0058,-0.0067,-0.0075,-0.0084,-0.0092,-0.0099,-0.0106,-0.0112,-0.0116,-0.0120,-0.0122,-0.0123,-0.0123,-0.0121,-0.0118,-0.0113,-0.0108,-0.0101,-0.0095,-0.0089,-0.0084,-0.0080,-0.0077,-0.0075,-0.0074,
                                                                0.0237,0.0177,0.0129,0.0091,0.0062,0.0040,0.0025,0.0015 ,0.0009 ,0.0004 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0006,-0.0008,-0.0010,-0.0012,-0.0014,-0.0017,-0.0021,-0.0025,-0.0030,-0.0035,-0.0042,-0.0049,-0.0057,-0.0065,-0.0074,-0.0082,-0.0091,-0.0098,-0.0105,-0.0110,-0.0115,-0.0118,-0.0120,-0.0121,-0.0120,-0.0117,-0.0113,-0.0108,-0.0102,-0.0096,-0.0090,-0.0085,-0.0081,-0.0078,-0.0075,-0.0074,
                                                                0.0253,0.0190,0.0139,0.0099,0.0068,0.0046,0.0030,0.0019 ,0.0011 ,0.0006 ,0.0003 ,0.0001 ,-0.0001,-0.0003,-0.0005,-0.0006,-0.0008,-0.0010,-0.0012,-0.0014,-0.0017,-0.0020,-0.0024,-0.0029,-0.0035,-0.0041,-0.0048,-0.0056,-0.0064,-0.0073,-0.0081,-0.0089,-0.0097,-0.0103,-0.0109,-0.0113,-0.0116,-0.0118,-0.0118,-0.0116,-0.0113,-0.0108,-0.0102,-0.0096,-0.0090,-0.0085,-0.0081,-0.0078,-0.0076,-0.0074,
                                                                0.0269,0.0204,0.0150,0.0108,0.0075,0.0051,0.0034,0.0022 ,0.0014 ,0.0009 ,0.0005 ,0.0002 ,-0.0000,-0.0002,-0.0003,-0.0005,-0.0007,-0.0008,-0.0010,-0.0012,-0.0014,-0.0017,-0.0020,-0.0024,-0.0029,-0.0034,-0.0040,-0.0047,-0.0055,-0.0063,-0.0071,-0.0080,-0.0088,-0.0095,-0.0102,-0.0107,-0.0111,-0.0114,-0.0115,-0.0114,-0.0112,-0.0108,-0.0102,-0.0097,-0.0091,-0.0086,-0.0081,-0.0078,-0.0076,-0.0074,
                                                                0.0287,0.0219,0.0163,0.0118,0.0083,0.0058,0.0039,0.0026 ,0.0017 ,0.0011 ,0.0007 ,0.0004 ,0.0001 ,-0.0001,-0.0002,-0.0004,-0.0005,-0.0007,-0.0008,-0.0010,-0.0012,-0.0014,-0.0017,-0.0020,-0.0024,-0.0028,-0.0033,-0.0039,-0.0046,-0.0054,-0.0062,-0.0070,-0.0078,-0.0086,-0.0094,-0.0100,-0.0105,-0.0109,-0.0111,-0.0111,-0.0110,-0.0107,-0.0102,-0.0097,-0.0091,-0.0086,-0.0082,-0.0078,-0.0076,-0.0074,
                                                                0.0306,0.0235,0.0176,0.0129,0.0092,0.0064,0.0044,0.0030 ,0.0020 ,0.0014 ,0.0009 ,0.0005 ,0.0003 ,0.0000 ,-0.0001,-0.0003,-0.0004,-0.0006,-0.0007,-0.0009,-0.0010,-0.0012,-0.0014,-0.0017,-0.0020,-0.0023,-0.0027,-0.0032,-0.0038,-0.0045,-0.0052,-0.0060,-0.0069,-0.0077,-0.0085,-0.0092,-0.0098,-0.0103,-0.0106,-0.0108,-0.0107,-0.0105,-0.0101,-0.0096,-0.0091,-0.0086,-0.0082,-0.0078,-0.0076,-0.0074,
                                                                0.0325,0.0252,0.0190,0.0141,0.0102,0.0072,0.0050,0.0035 ,0.0024 ,0.0016 ,0.0011 ,0.0007 ,0.0004 ,0.0002 ,-0.0000,-0.0002,-0.0004,-0.0005,-0.0006,-0.0008,-0.0009,-0.0010,-0.0012,-0.0014,-0.0016,-0.0019,-0.0023,-0.0027,-0.0032,-0.0037,-0.0044,-0.0051,-0.0059,-0.0067,-0.0075,-0.0083,-0.0090,-0.0096,-0.0100,-0.0103,-0.0104,-0.0102,-0.0099,-0.0095,-0.0090,-0.0086,-0.0082,-0.0078,-0.0076,-0.0074,
                                                                0.0346,0.0270,0.0206,0.0154,0.0112,0.0081,0.0057,0.0040 ,0.0028 ,0.0020 ,0.0013 ,0.0009 ,0.0006 ,0.0003 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0005,-0.0007,-0.0008,-0.0009,-0.0010,-0.0012,-0.0014,-0.0016,-0.0019,-0.0022,-0.0026,-0.0031,-0.0037,-0.0043,-0.0050,-0.0058,-0.0066,-0.0074,-0.0081,-0.0088,-0.0093,-0.0097,-0.0099,-0.0099,-0.0097,-0.0093,-0.0089,-0.0085,-0.0081,-0.0078,-0.0075,-0.0074,
                                                                0.0367,0.0289,0.0223,0.0168,0.0124,0.0090,0.0065,0.0046 ,0.0033 ,0.0023 ,0.0016 ,0.0011 ,0.0007 ,0.0004 ,0.0002 ,-0.0000,-0.0002,-0.0004,-0.0005,-0.0006,-0.0007,-0.0008,-0.0009,-0.0011,-0.0012,-0.0014,-0.0016,-0.0019,-0.0022,-0.0026,-0.0030,-0.0036,-0.0042,-0.0049,-0.0056,-0.0064,-0.0072,-0.0079,-0.0085,-0.0090,-0.0093,-0.0094,-0.0093,-0.0091,-0.0087,-0.0084,-0.0080,-0.0077,-0.0075,-0.0073,
                                                                0.0389,0.0309,0.0240,0.0183,0.0137,0.0101,0.0073,0.0053 ,0.0038 ,0.0027 ,0.0019 ,0.0013 ,0.0009 ,0.0006 ,0.0003 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0005,-0.0006,-0.0007,-0.0008,-0.0009,-0.0011,-0.0012,-0.0014,-0.0016,-0.0018,-0.0021,-0.0025,-0.0029,-0.0035,-0.0041,-0.0047,-0.0055,-0.0062,-0.0069,-0.0076,-0.0082,-0.0086,-0.0088,-0.0088,-0.0087,-0.0085,-0.0081,-0.0078,-0.0076,-0.0074,-0.0072,
                                                                0.0410,0.0329,0.0259,0.0199,0.0151,0.0112,0.0083,0.0060 ,0.0044 ,0.0032 ,0.0023 ,0.0016 ,0.0011 ,0.0007 ,0.0004 ,0.0002 ,-0.0001,-0.0002,-0.0004,-0.0005,-0.0006,-0.0007,-0.0007,-0.0008,-0.0009,-0.0011,-0.0012,-0.0013,-0.0015,-0.0018,-0.0021,-0.0024,-0.0029,-0.0034,-0.0039,-0.0046,-0.0053,-0.0060,-0.0067,-0.0073,-0.0078,-0.0081,-0.0083,-0.0082,-0.0081,-0.0079,-0.0076,-0.0074,-0.0072,-0.0071,
                                                                0.0432,0.0349,0.0277,0.0216,0.0166,0.0125,0.0093,0.0069 ,0.0051 ,0.0037 ,0.0027 ,0.0019 ,0.0014 ,0.0009 ,0.0005 ,0.0003 ,0.0000 ,-0.0001,-0.0003,-0.0004,-0.0005,-0.0006,-0.0007,-0.0008,-0.0009,-0.0009,-0.0010,-0.0012,-0.0013,-0.0015,-0.0017,-0.0020,-0.0024,-0.0028,-0.0033,-0.0038,-0.0044,-0.0051,-0.0057,-0.0063,-0.0069,-0.0073,-0.0076,-0.0077,-0.0076,-0.0075,-0.0073,-0.0071,-0.0070,-0.0069,
                                                                0.0453,0.0370,0.0296,0.0234,0.0181,0.0139,0.0105,0.0078 ,0.0058 ,0.0043 ,0.0032 ,0.0023 ,0.0016 ,0.0011 ,0.0007 ,0.0004 ,0.0001 ,-0.0001,-0.0002,-0.0004,-0.0005,-0.0006,-0.0006,-0.0007,-0.0008,-0.0009,-0.0009,-0.0010,-0.0012,-0.0013,-0.0015,-0.0017,-0.0020,-0.0023,-0.0027,-0.0031,-0.0037,-0.0042,-0.0048,-0.0054,-0.0060,-0.0064,-0.0068,-0.0070,-0.0070,-0.0070,-0.0069,-0.0068,-0.0067,-0.0067,
                                                                0.0473,0.0389,0.0315,0.0252,0.0198,0.0153,0.0117,0.0089 ,0.0067 ,0.0050 ,0.0037 ,0.0027 ,0.0020 ,0.0014 ,0.0009 ,0.0005 ,0.0002 ,0.0000 ,-0.0002,-0.0003,-0.0004,-0.0005,-0.0006,-0.0007,-0.0007,-0.0008,-0.0009,-0.0009,-0.0010,-0.0012,-0.0013,-0.0015,-0.0017,-0.0019,-0.0022,-0.0026,-0.0030,-0.0035,-0.0040,-0.0046,-0.0051,-0.0056,-0.0059,-0.0062,-0.0063,-0.0063,-0.0063,-0.0063,-0.0063,-0.0064,
                                                                0.0493,0.0409,0.0334,0.0269,0.0214,0.0168,0.0131,0.0100 ,0.0077 ,0.0058 ,0.0043 ,0.0032 ,0.0023 ,0.0017 ,0.0011 ,0.0007 ,0.0004 ,0.0001 ,-0.0001,-0.0002,-0.0004,-0.0005,-0.0005,-0.0006,-0.0007,-0.0007,-0.0008,-0.0009,-0.0009,-0.0010,-0.0011,-0.0013,-0.0014,-0.0016,-0.0019,-0.0021,-0.0025,-0.0029,-0.0033,-0.0038,-0.0042,-0.0047,-0.0050,-0.0053,-0.0055,-0.0056,-0.0057,-0.0058,-0.0058,-0.0059,
                                                                0.0512,0.0428,0.0353,0.0287,0.0231,0.0184,0.0145,0.0113 ,0.0087 ,0.0067 ,0.0051 ,0.0038 ,0.0028 ,0.0020 ,0.0014 ,0.0009 ,0.0005 ,0.0002 ,0.0000 ,-0.0002,-0.0003,-0.0004,-0.0005,-0.0006,-0.0006,-0.0007,-0.0008,-0.0008,-0.0009,-0.0009,-0.0010,-0.0011,-0.0012,-0.0014,-0.0016,-0.0018,-0.0021,-0.0024,-0.0027,-0.0031,-0.0035,-0.0039,-0.0042,-0.0045,-0.0047,-0.0048,-0.0050,-0.0051,-0.0053,-0.0054,
                                                                0.0531,0.0446,0.0370,0.0304,0.0248,0.0200,0.0159,0.0126 ,0.0099 ,0.0077 ,0.0059 ,0.0044 ,0.0033 ,0.0024 ,0.0017 ,0.0011 ,0.0007 ,0.0004 ,0.0001 ,-0.0001,-0.0002,-0.0004,-0.0005,-0.0005,-0.0006,-0.0007,-0.0007,-0.0008,-0.0008,-0.0009,-0.0009,-0.0010,-0.0011,-0.0012,-0.0013,-0.0015,-0.0017,-0.0019,-0.0022,-0.0025,-0.0028,-0.0031,-0.0034,-0.0037,-0.0039,-0.0040,-0.0042,-0.0043,-0.0046,-0.0048,
                                                                0.0548,0.0463,0.0387,0.0321,0.0264,0.0215,0.0174,0.0140 ,0.0111 ,0.0087 ,0.0068 ,0.0052 ,0.0039 ,0.0029 ,0.0020 ,0.0014 ,0.0009 ,0.0005 ,0.0002 ,-0.0000,-0.0002,-0.0003,-0.0004,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0008,-0.0009,-0.0009,-0.0010,-0.0011,-0.0012,-0.0013,-0.0014,-0.0016,-0.0018,-0.0020,-0.0023,-0.0025,-0.0027,-0.0029,-0.0030,-0.0032,-0.0033,-0.0035,-0.0038,-0.0041,
                                                                0.0566,0.0480,0.0404,0.0337,0.0280,0.0231,0.0189,0.0154 ,0.0124 ,0.0099 ,0.0078 ,0.0060 ,0.0046 ,0.0034 ,0.0025 ,0.0017 ,0.0012 ,0.0007 ,0.0004 ,0.0001 ,-0.0001,-0.0002,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0008,-0.0009,-0.0009,-0.0010,-0.0010,-0.0011,-0.0012,-0.0014,-0.0015,-0.0017,-0.0018,-0.0020,-0.0021,-0.0022,-0.0023,-0.0024,-0.0025,-0.0027,-0.0029,-0.0033,
                                                                0.0583,0.0496,0.0419,0.0352,0.0295,0.0246,0.0204,0.0168 ,0.0138 ,0.0111 ,0.0089 ,0.0070 ,0.0054 ,0.0040 ,0.0030 ,0.0021 ,0.0015 ,0.0009 ,0.0005 ,0.0002 ,0.0000 ,-0.0002,-0.0003,-0.0004,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0007,-0.0008,-0.0008,-0.0008,-0.0009,-0.0010,-0.0010,-0.0011,-0.0012,-0.0013,-0.0014,-0.0015,-0.0015,-0.0016,-0.0016,-0.0017,-0.0017,-0.0017,-0.0018,-0.0021,-0.0024,
                                                                0.0599,0.0512,0.0435,0.0367,0.0309,0.0260,0.0218,0.0182 ,0.0151 ,0.0124 ,0.0100 ,0.0080 ,0.0062 ,0.0048 ,0.0036 ,0.0026 ,0.0018 ,0.0012 ,0.0008 ,0.0004 ,0.0001 ,-0.0001,-0.0002,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0007,-0.0008,-0.0008,-0.0008,-0.0009,-0.0009,-0.0010,-0.0010,-0.0011,-0.0011,-0.0012,-0.0012,-0.0012,-0.0012,-0.0011,-0.0011,-0.0010,-0.0011,-0.0012,-0.0015,
                                                                0.0616,0.0528,0.0449,0.0381,0.0323,0.0273,0.0231,0.0196 ,0.0164 ,0.0137 ,0.0112 ,0.0091 ,0.0072 ,0.0056 ,0.0042 ,0.0031 ,0.0023 ,0.0016 ,0.0010 ,0.0006 ,0.0003 ,0.0000 ,-0.0002,-0.0003,-0.0004,-0.0005,-0.0006,-0.0006,-0.0006,-0.0007,-0.0007,-0.0007,-0.0008,-0.0008,-0.0008,-0.0009,-0.0009,-0.0009,-0.0010,-0.0010,-0.0010,-0.0010,-0.0009,-0.0008,-0.0007,-0.0005,-0.0004,-0.0004,-0.0004,-0.0007,
                                                                0.0632,0.0543,0.0464,0.0395,0.0336,0.0286,0.0244,0.0209 ,0.0177 ,0.0150 ,0.0125 ,0.0102 ,0.0082 ,0.0065 ,0.0050 ,0.0038 ,0.0028 ,0.0019 ,0.0013 ,0.0008 ,0.0004 ,0.0002 ,-0.0001,-0.0002,-0.0004,-0.0004,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0007,-0.0007,-0.0008,-0.0008,-0.0008,-0.0008,-0.0008,-0.0009,-0.0009,-0.0008,-0.0008,-0.0007,-0.0006,-0.0004,-0.0001,0.0001 ,0.0002 ,0.0003 ,0.0001 ,
                                                                0.0648,0.0558,0.0478,0.0408,0.0349,0.0298,0.0256,0.0221 ,0.0190 ,0.0162 ,0.0137 ,0.0114 ,0.0094 ,0.0075 ,0.0059 ,0.0045 ,0.0033 ,0.0024 ,0.0017 ,0.0011 ,0.0007 ,0.0003 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0008,-0.0008,-0.0008,-0.0008,-0.0008,-0.0007,-0.0007,-0.0005,-0.0004,-0.0001,0.0002 ,0.0005 ,0.0007 ,0.0009 ,0.0008 ,
                                                                0.0664,0.0573,0.0492,0.0421,0.0361,0.0310,0.0268,0.0232 ,0.0201 ,0.0174 ,0.0149 ,0.0126 ,0.0105 ,0.0086 ,0.0068 ,0.0053 ,0.0040 ,0.0030 ,0.0021 ,0.0014 ,0.0009 ,0.0005 ,0.0002 ,-0.0000,-0.0002,-0.0003,-0.0004,-0.0005,-0.0006,-0.0006,-0.0006,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0006,-0.0004,-0.0002,0.0000 ,0.0004 ,0.0007 ,0.0011 ,0.0013 ,0.0014 ,
                                                                0.0679,0.0588,0.0506,0.0434,0.0372,0.0321,0.0279,0.0243 ,0.0212 ,0.0186 ,0.0161 ,0.0138 ,0.0117 ,0.0097 ,0.0078 ,0.0062 ,0.0048 ,0.0036 ,0.0026 ,0.0018 ,0.0012 ,0.0007 ,0.0004 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0006,-0.0005,-0.0004,-0.0001,0.0002 ,0.0005 ,0.0009 ,0.0013 ,0.0017 ,0.0019 ,
                                                                0.0695,0.0603,0.0519,0.0446,0.0384,0.0332,0.0289,0.0253 ,0.0223 ,0.0196 ,0.0172 ,0.0150 ,0.0128 ,0.0108 ,0.0089 ,0.0072 ,0.0056 ,0.0043 ,0.0032 ,0.0023 ,0.0016 ,0.0010 ,0.0006 ,0.0003 ,0.0000 ,-0.0002,-0.0003,-0.0004,-0.0005,-0.0006,-0.0006,-0.0006,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0006,-0.0006,-0.0005,-0.0003,-0.0001,0.0002 ,0.0006 ,0.0010 ,0.0015 ,0.0019 ,0.0023 ,
                                                                0.0711,0.0618,0.0533,0.0459,0.0395,0.0342,0.0298,0.0262 ,0.0232 ,0.0206 ,0.0182 ,0.0161 ,0.0140 ,0.0120 ,0.0100 ,0.0082 ,0.0066 ,0.0051 ,0.0039 ,0.0029 ,0.0020 ,0.0014 ,0.0009 ,0.0005 ,0.0002 ,-0.0001,-0.0002,-0.0004,-0.0004,-0.0005,-0.0006,-0.0006,-0.0006,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0006,-0.0006,-0.0005,-0.0004,-0.0003,-0.0001,0.0002 ,0.0006 ,0.0011 ,0.0016 ,0.0021 ,0.0025 ,
                                                                0.0727,0.0632,0.0547,0.0471,0.0406,0.0352,0.0308,0.0271 ,0.0241 ,0.0215 ,0.0192 ,0.0171 ,0.0150 ,0.0131 ,0.0112 ,0.0093 ,0.0076 ,0.0060 ,0.0047 ,0.0035 ,0.0026 ,0.0018 ,0.0012 ,0.0007 ,0.0004 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0006,-0.0007,-0.0007,-0.0007,-0.0007,-0.0006,-0.0006,-0.0005,-0.0004,-0.0003,-0.0001,0.0003 ,0.0006 ,0.0011 ,0.0017 ,0.0022 ,0.0027 ,
                                                                0.0743,0.0647,0.0560,0.0484,0.0418,0.0362,0.0317,0.0280 ,0.0249 ,0.0223 ,0.0201 ,0.0180 ,0.0160 ,0.0141 ,0.0123 ,0.0104 ,0.0087 ,0.0070 ,0.0056 ,0.0043 ,0.0032 ,0.0023 ,0.0016 ,0.0010 ,0.0006 ,0.0003 ,0.0000 ,-0.0002,-0.0003,-0.0004,-0.0005,-0.0006,-0.0006,-0.0006,-0.0006,-0.0007,-0.0007,-0.0006,-0.0006,-0.0006,-0.0005,-0.0004,-0.0003,-0.0001,0.0002 ,0.0006 ,0.0011 ,0.0017 ,0.0023 ,0.0028 ,
                                                                0.0759,0.0663,0.0574,0.0496,0.0429,0.0372,0.0326,0.0288 ,0.0257 ,0.0231 ,0.0209 ,0.0188 ,0.0170 ,0.0151 ,0.0133 ,0.0115 ,0.0097 ,0.0081 ,0.0065 ,0.0051 ,0.0039 ,0.0029 ,0.0020 ,0.0014 ,0.0009 ,0.0005 ,0.0002 ,-0.0001,-0.0002,-0.0004,-0.0004,-0.0005,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0005,-0.0004,-0.0003,-0.0001,0.0002 ,0.0006 ,0.0011 ,0.0016 ,0.0022 ,0.0028 ,
                                                                0.0775,0.0678,0.0588,0.0509,0.0440,0.0382,0.0335,0.0296 ,0.0264 ,0.0238 ,0.0216 ,0.0196 ,0.0178 ,0.0160 ,0.0143 ,0.0126 ,0.0108 ,0.0091 ,0.0075 ,0.0060 ,0.0047 ,0.0035 ,0.0026 ,0.0018 ,0.0012 ,0.0007 ,0.0004 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0005,-0.0004,-0.0003,-0.0001,0.0002 ,0.0005 ,0.0010 ,0.0016 ,0.0022 ,0.0028 ,
                                                                0.0791,0.0693,0.0603,0.0522,0.0451,0.0392,0.0343,0.0304 ,0.0271 ,0.0245 ,0.0222 ,0.0203 ,0.0185 ,0.0169 ,0.0152 ,0.0136 ,0.0119 ,0.0102 ,0.0086 ,0.0070 ,0.0056 ,0.0043 ,0.0032 ,0.0023 ,0.0016 ,0.0010 ,0.0006 ,0.0003 ,0.0000 ,-0.0002,-0.0003,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0005,-0.0004,-0.0003,-0.0001,0.0001 ,0.0005 ,0.0010 ,0.0015 ,0.0021 ,0.0028 ,
                                                                0.0807,0.0708,0.0617,0.0534,0.0463,0.0402,0.0352,0.0311 ,0.0278 ,0.0251 ,0.0229 ,0.0209 ,0.0192 ,0.0176 ,0.0160 ,0.0145 ,0.0129 ,0.0113 ,0.0096 ,0.0080 ,0.0065 ,0.0051 ,0.0039 ,0.0029 ,0.0021 ,0.0014 ,0.0009 ,0.0005 ,0.0002 ,-0.0000,-0.0002,-0.0003,-0.0004,-0.0005,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0005,-0.0005,-0.0003,-0.0002,0.0001 ,0.0004 ,0.0009 ,0.0014 ,0.0021 ,0.0027 ,
                                                                0.0823,0.0724,0.0631,0.0548,0.0474,0.0412,0.0361,0.0319 ,0.0285 ,0.0257 ,0.0234 ,0.0215 ,0.0198 ,0.0183 ,0.0168 ,0.0153 ,0.0138 ,0.0123 ,0.0107 ,0.0091 ,0.0075 ,0.0061 ,0.0048 ,0.0036 ,0.0027 ,0.0019 ,0.0013 ,0.0008 ,0.0004 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0006,-0.0006,-0.0006,-0.0005,-0.0005,-0.0004,-0.0002,0.0000 ,0.0004 ,0.0008 ,0.0013 ,0.0020 ,0.0026 ,
                                                                0.0839,0.0740,0.0646,0.0561,0.0486,0.0423,0.0370,0.0327 ,0.0291 ,0.0263 ,0.0240 ,0.0220 ,0.0204 ,0.0188 ,0.0174 ,0.0161 ,0.0146 ,0.0132 ,0.0117 ,0.0101 ,0.0086 ,0.0071 ,0.0057 ,0.0044 ,0.0033 ,0.0024 ,0.0017 ,0.0011 ,0.0007 ,0.0003 ,0.0000 ,-0.0002,-0.0003,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0006,-0.0006,-0.0005,-0.0005,-0.0004,-0.0002,-0.0000,0.0003 ,0.0007 ,0.0012 ,0.0019 ,0.0025 ,
                                                                0.0856,0.0755,0.0661,0.0574,0.0498,0.0433,0.0379,0.0334 ,0.0298 ,0.0269 ,0.0245 ,0.0225 ,0.0208 ,0.0194 ,0.0180 ,0.0167 ,0.0154 ,0.0141 ,0.0126 ,0.0111 ,0.0096 ,0.0081 ,0.0066 ,0.0053 ,0.0041 ,0.0030 ,0.0022 ,0.0015 ,0.0010 ,0.0005 ,0.0002 ,-0.0000,-0.0002,-0.0003,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0006,-0.0005,-0.0005,-0.0004,-0.0003,-0.0000,0.0002 ,0.0006 ,0.0011 ,0.0017 ,0.0024 ,
                                                                0.0872,0.0771,0.0676,0.0588,0.0510,0.0444,0.0388,0.0342 ,0.0305 ,0.0275 ,0.0250 ,0.0230 ,0.0213 ,0.0198 ,0.0185 ,0.0173 ,0.0161 ,0.0148 ,0.0135 ,0.0121 ,0.0106 ,0.0091 ,0.0076 ,0.0062 ,0.0049 ,0.0038 ,0.0028 ,0.0020 ,0.0013 ,0.0008 ,0.0004 ,0.0001 ,-0.0001,-0.0002,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0005,-0.0005,-0.0004,-0.0003,-0.0001,0.0002 ,0.0006 ,0.0010 ,0.0016 ,0.0023 ,
                                                                0.0888,0.0787,0.0691,0.0602,0.0523,0.0455,0.0397,0.0350 ,0.0311 ,0.0280 ,0.0255 ,0.0234 ,0.0217 ,0.0203 ,0.0190 ,0.0178 ,0.0167 ,0.0155 ,0.0143 ,0.0130 ,0.0116 ,0.0102 ,0.0087 ,0.0072 ,0.0058 ,0.0046 ,0.0035 ,0.0025 ,0.0018 ,0.0012 ,0.0007 ,0.0004 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0005,-0.0005,-0.0005,-0.0005,-0.0005,-0.0005,-0.0004,-0.0003,-0.0001,0.0001 ,0.0005 ,0.0009 ,0.0015 ,0.0022 ,
                                                                0.0904,0.0803,0.0706,0.0616,0.0535,0.0466,0.0407,0.0358 ,0.0318 ,0.0286 ,0.0260 ,0.0239 ,0.0221 ,0.0207 ,0.0194 ,0.0183 ,0.0172 ,0.0161 ,0.0150 ,0.0138 ,0.0125 ,0.0111 ,0.0097 ,0.0082 ,0.0068 ,0.0054 ,0.0042 ,0.0032 ,0.0023 ,0.0016 ,0.0010 ,0.0006 ,0.0003 ,0.0000 ,-0.0002,-0.0003,-0.0004,-0.0005,-0.0005,-0.0005,-0.0005,-0.0005,-0.0004,-0.0003,-0.0002,0.0001 ,0.0004 ,0.0008 ,0.0014 ,0.0020 ,
                                                                0.0921,0.0819,0.0721,0.0630,0.0548,0.0477,0.0416,0.0366 ,0.0325 ,0.0292 ,0.0265 ,0.0243 ,0.0225 ,0.0210 ,0.0198 ,0.0186 ,0.0176 ,0.0166 ,0.0156 ,0.0145 ,0.0134 ,0.0121 ,0.0107 ,0.0093 ,0.0078 ,0.0064 ,0.0051 ,0.0039 ,0.0029 ,0.0021 ,0.0014 ,0.0009 ,0.0005 ,0.0002 ,-0.0000,-0.0002,-0.0003,-0.0004,-0.0005,-0.0005,-0.0005,-0.0005,-0.0004,-0.0003,-0.0002,0.0000 ,0.0003 ,0.0008 ,0.0013 ,0.0019 ,
                                                                0.0937,0.0835,0.0737,0.0644,0.0561,0.0488,0.0426,0.0375 ,0.0332 ,0.0298 ,0.0270 ,0.0247 ,0.0229 ,0.0214 ,0.0201 ,0.0190 ,0.0180 ,0.0171 ,0.0161 ,0.0152 ,0.0141 ,0.0129 ,0.0116 ,0.0103 ,0.0088 ,0.0074 ,0.0060 ,0.0048 ,0.0036 ,0.0027 ,0.0019 ,0.0013 ,0.0008 ,0.0004 ,0.0001 ,-0.0001,-0.0002,-0.0004,-0.0004,-0.0005,-0.0005,-0.0005,-0.0004,-0.0004,-0.0002,-0.0000,0.0003 ,0.0007 ,0.0012 ,0.0018 ,
                                                                0.0953,0.0851,0.0752,0.0659,0.0574,0.0500,0.0436,0.0383 ,0.0339 ,0.0304 ,0.0275 ,0.0251 ,0.0232 ,0.0217 ,0.0204 ,0.0193 ,0.0184 ,0.0175 ,0.0166 ,0.0157 ,0.0148 ,0.0137 ,0.0125 ,0.0112 ,0.0099 ,0.0084 ,0.0070 ,0.0057 ,0.0044 ,0.0034 ,0.0025 ,0.0017 ,0.0011 ,0.0007 ,0.0003 ,0.0001 ,-0.0001,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0004,-0.0004,-0.0002,-0.0001,0.0002 ,0.0006 ,0.0011 ,0.0016 ,
                                                                0.0969,0.0868,0.0768,0.0674,0.0588,0.0512,0.0447,0.0392 ,0.0347 ,0.0310 ,0.0280 ,0.0256 ,0.0236 ,0.0220 ,0.0207 ,0.0196 ,0.0187 ,0.0178 ,0.0170 ,0.0162 ,0.0153 ,0.0144 ,0.0133 ,0.0121 ,0.0108 ,0.0094 ,0.0080 ,0.0066 ,0.0053 ,0.0041 ,0.0031 ,0.0022 ,0.0016 ,0.0010 ,0.0006 ,0.0003 ,0.0000 ,-0.0002,-0.0003,-0.0004,-0.0004,-0.0004,-0.0004,-0.0004,-0.0003,-0.0001,0.0001 ,0.0005 ,0.0010 ,0.0015};
        
        return arrayCurvatureMin[y * arrayWidth + x];
    }
    case 9:
    {
        static const std::vector<double> arrayCurvatureMin   = {0.0200,0.0153,0.0112,0.0079,0.0051,0.0030,0.0013,0.0001,-0.0007,-0.0012,-0.0015,-0.0016,-0.0016,-0.0015,-0.0014,-0.0013,-0.0012,-0.0012,-0.0012,-0.0012,-0.0013,-0.0015,-0.0017,-0.0020,-0.0023,-0.0027,-0.0033,-0.0038,-0.0045,-0.0052,-0.0060,-0.0067,-0.0074,-0.0080,-0.0086,-0.0089,-0.0091,-0.0091,-0.0088,-0.0084,-0.0078,-0.0070,-0.0061,-0.0052,-0.0042,-0.0032,-0.0022,-0.0012,-0.0003,0.0005,
                                                                0.0211,0.0163,0.0122,0.0088,0.0060,0.0038,0.0021,0.0008,-0.0000,-0.0006,-0.0009,-0.0011,-0.0011,-0.0011,-0.0011,-0.0010,-0.0010,-0.0010,-0.0010,-0.0011,-0.0012,-0.0014,-0.0016,-0.0019,-0.0023,-0.0027,-0.0032,-0.0038,-0.0045,-0.0052,-0.0059,-0.0066,-0.0073,-0.0079,-0.0084,-0.0087,-0.0088,-0.0087,-0.0085,-0.0080,-0.0074,-0.0066,-0.0057,-0.0048,-0.0038,-0.0028,-0.0019,-0.0009,-0.0001,0.0007,
                                                                0.0223,0.0175,0.0133,0.0098,0.0069,0.0046,0.0029,0.0016,0.0006 ,0.0000 ,-0.0004,-0.0006,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0007,-0.0008,-0.0009,-0.0011,-0.0012,-0.0015,-0.0018,-0.0022,-0.0026,-0.0032,-0.0037,-0.0044,-0.0051,-0.0058,-0.0065,-0.0071,-0.0077,-0.0081,-0.0084,-0.0085,-0.0084,-0.0081,-0.0076,-0.0070,-0.0062,-0.0053,-0.0044,-0.0035,-0.0025,-0.0016,-0.0007,0.0001 ,0.0009,
                                                                0.0236,0.0186,0.0143,0.0107,0.0078,0.0055,0.0037,0.0023,0.0013 ,0.0007 ,0.0002 ,-0.0001,-0.0002,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0008,-0.0009,-0.0011,-0.0014,-0.0017,-0.0021,-0.0026,-0.0031,-0.0037,-0.0043,-0.0050,-0.0057,-0.0064,-0.0070,-0.0075,-0.0079,-0.0081,-0.0082,-0.0080,-0.0077,-0.0072,-0.0066,-0.0058,-0.0050,-0.0041,-0.0031,-0.0022,-0.0013,-0.0004,0.0004 ,0.0011,
                                                                0.0249,0.0198,0.0154,0.0118,0.0088,0.0064,0.0045,0.0031,0.0020 ,0.0013 ,0.0008 ,0.0005 ,0.0002 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0003,-0.0005,-0.0006,-0.0008,-0.0010,-0.0013,-0.0016,-0.0020,-0.0025,-0.0030,-0.0036,-0.0042,-0.0049,-0.0056,-0.0062,-0.0068,-0.0073,-0.0076,-0.0078,-0.0078,-0.0076,-0.0073,-0.0068,-0.0062,-0.0054,-0.0046,-0.0037,-0.0028,-0.0019,-0.0011,-0.0002,0.0006 ,0.0013,
                                                                0.0262,0.0210,0.0166,0.0128,0.0098,0.0073,0.0053,0.0039,0.0028 ,0.0019 ,0.0014 ,0.0010 ,0.0007 ,0.0005 ,0.0003 ,0.0001 ,0.0000 ,-0.0001,-0.0003,-0.0005,-0.0007,-0.0009,-0.0012,-0.0015,-0.0019,-0.0024,-0.0029,-0.0035,-0.0041,-0.0048,-0.0054,-0.0060,-0.0066,-0.0070,-0.0073,-0.0075,-0.0074,-0.0073,-0.0069,-0.0064,-0.0058,-0.0050,-0.0042,-0.0034,-0.0025,-0.0017,-0.0008,-0.0000,0.0007 ,0.0014,
                                                                0.0275,0.0223,0.0178,0.0139,0.0108,0.0082,0.0062,0.0046,0.0035 ,0.0026 ,0.0020 ,0.0015 ,0.0011 ,0.0008 ,0.0006 ,0.0004 ,0.0002 ,0.0001 ,-0.0001,-0.0003,-0.0005,-0.0008,-0.0011,-0.0014,-0.0018,-0.0023,-0.0028,-0.0034,-0.0040,-0.0046,-0.0052,-0.0058,-0.0063,-0.0067,-0.0070,-0.0071,-0.0071,-0.0069,-0.0065,-0.0060,-0.0054,-0.0047,-0.0039,-0.0031,-0.0023,-0.0014,-0.0006,0.0002 ,0.0009 ,0.0015,
                                                                0.0289,0.0236,0.0190,0.0151,0.0118,0.0092,0.0071,0.0055,0.0042 ,0.0033 ,0.0025 ,0.0020 ,0.0015 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0003 ,0.0001 ,-0.0002,-0.0004,-0.0007,-0.0010,-0.0013,-0.0017,-0.0022,-0.0027,-0.0033,-0.0038,-0.0045,-0.0050,-0.0056,-0.0061,-0.0064,-0.0067,-0.0068,-0.0067,-0.0065,-0.0061,-0.0056,-0.0050,-0.0043,-0.0036,-0.0028,-0.0020,-0.0012,-0.0004,0.0003 ,0.0010 ,0.0017,
                                                                0.0303,0.0249,0.0202,0.0162,0.0129,0.0102,0.0080,0.0063,0.0049 ,0.0039 ,0.0031 ,0.0025 ,0.0020 ,0.0016 ,0.0012 ,0.0010 ,0.0007 ,0.0005 ,0.0002 ,0.0000 ,-0.0003,-0.0005,-0.0008,-0.0012,-0.0016,-0.0021,-0.0026,-0.0031,-0.0037,-0.0043,-0.0048,-0.0054,-0.0058,-0.0061,-0.0063,-0.0064,-0.0063,-0.0061,-0.0057,-0.0052,-0.0046,-0.0040,-0.0033,-0.0025,-0.0018,-0.0010,-0.0003,0.0005 ,0.0011 ,0.0017,
                                                                0.0318,0.0263,0.0215,0.0174,0.0140,0.0112,0.0089,0.0071,0.0057 ,0.0046 ,0.0037 ,0.0030 ,0.0024 ,0.0019 ,0.0016 ,0.0012 ,0.0009 ,0.0007 ,0.0004 ,0.0002 ,-0.0001,-0.0004,-0.0007,-0.0011,-0.0015,-0.0019,-0.0024,-0.0029,-0.0035,-0.0041,-0.0046,-0.0051,-0.0055,-0.0058,-0.0060,-0.0060,-0.0059,-0.0057,-0.0053,-0.0048,-0.0043,-0.0036,-0.0030,-0.0023,-0.0015,-0.0008,-0.0001,0.0006 ,0.0012 ,0.0018,
                                                                0.0333,0.0277,0.0229,0.0187,0.0152,0.0122,0.0099,0.0080,0.0064 ,0.0052 ,0.0042 ,0.0035 ,0.0028 ,0.0023 ,0.0019 ,0.0015 ,0.0012 ,0.0009 ,0.0006 ,0.0003 ,0.0001 ,-0.0002,-0.0006,-0.0009,-0.0013,-0.0018,-0.0022,-0.0028,-0.0033,-0.0038,-0.0044,-0.0048,-0.0052,-0.0055,-0.0056,-0.0057,-0.0055,-0.0053,-0.0049,-0.0045,-0.0039,-0.0033,-0.0027,-0.0020,-0.0013,-0.0006,0.0000 ,0.0007 ,0.0013 ,0.0019,
                                                                0.0348,0.0292,0.0242,0.0199,0.0163,0.0133,0.0108,0.0088,0.0072 ,0.0059 ,0.0048 ,0.0039 ,0.0032 ,0.0027 ,0.0022 ,0.0018 ,0.0014 ,0.0011 ,0.0008 ,0.0005 ,0.0002 ,-0.0001,-0.0004,-0.0008,-0.0011,-0.0016,-0.0021,-0.0026,-0.0031,-0.0036,-0.0041,-0.0045,-0.0049,-0.0051,-0.0053,-0.0053,-0.0052,-0.0049,-0.0046,-0.0041,-0.0036,-0.0030,-0.0024,-0.0018,-0.0011,-0.0005,0.0002 ,0.0008 ,0.0014 ,0.0019,
                                                                0.0364,0.0307,0.0256,0.0212,0.0175,0.0144,0.0118,0.0097,0.0080 ,0.0065 ,0.0054 ,0.0044 ,0.0037 ,0.0030 ,0.0025 ,0.0020 ,0.0016 ,0.0013 ,0.0010 ,0.0007 ,0.0004 ,0.0001 ,-0.0002,-0.0006,-0.0010,-0.0014,-0.0019,-0.0023,-0.0028,-0.0033,-0.0038,-0.0042,-0.0046,-0.0048,-0.0049,-0.0049,-0.0048,-0.0045,-0.0042,-0.0038,-0.0033,-0.0028,-0.0022,-0.0016,-0.0010,-0.0004,0.0003 ,0.0009 ,0.0014 ,0.0020,
                                                                0.0380,0.0322,0.0270,0.0226,0.0187,0.0155,0.0128,0.0106,0.0087 ,0.0072 ,0.0060 ,0.0049 ,0.0041 ,0.0034 ,0.0028 ,0.0023 ,0.0019 ,0.0015 ,0.0012 ,0.0009 ,0.0006 ,0.0003 ,-0.0000,-0.0004,-0.0008,-0.0012,-0.0016,-0.0021,-0.0026,-0.0031,-0.0035,-0.0039,-0.0042,-0.0044,-0.0045,-0.0045,-0.0044,-0.0042,-0.0039,-0.0035,-0.0030,-0.0025,-0.0020,-0.0014,-0.0008,-0.0002,0.0004 ,0.0009 ,0.0015 ,0.0020,
                                                                0.0396,0.0337,0.0285,0.0239,0.0200,0.0166,0.0138,0.0114,0.0095 ,0.0079 ,0.0065 ,0.0054 ,0.0045 ,0.0037 ,0.0031 ,0.0026 ,0.0021 ,0.0017 ,0.0014 ,0.0011 ,0.0008 ,0.0005 ,0.0002 ,-0.0002,-0.0006,-0.0010,-0.0014,-0.0019,-0.0023,-0.0028,-0.0032,-0.0036,-0.0039,-0.0041,-0.0042,-0.0042,-0.0040,-0.0038,-0.0035,-0.0031,-0.0027,-0.0023,-0.0017,-0.0012,-0.0007,-0.0001,0.0004 ,0.0010 ,0.0015 ,0.0020,
                                                                0.0413,0.0353,0.0300,0.0253,0.0212,0.0177,0.0148,0.0123,0.0103 ,0.0085 ,0.0071 ,0.0059 ,0.0049 ,0.0041 ,0.0034 ,0.0029 ,0.0024 ,0.0020 ,0.0016 ,0.0013 ,0.0010 ,0.0007 ,0.0004 ,0.0000 ,-0.0003,-0.0007,-0.0011,-0.0016,-0.0020,-0.0025,-0.0029,-0.0032,-0.0035,-0.0037,-0.0038,-0.0038,-0.0037,-0.0035,-0.0032,-0.0029,-0.0025,-0.0020,-0.0016,-0.0011,-0.0005,-0.0000,0.0005 ,0.0010 ,0.0015 ,0.0020,
                                                                0.0430,0.0369,0.0315,0.0267,0.0225,0.0189,0.0158,0.0132,0.0110 ,0.0092 ,0.0077 ,0.0064 ,0.0053 ,0.0045 ,0.0038 ,0.0032 ,0.0027 ,0.0022 ,0.0019 ,0.0015 ,0.0012 ,0.0009 ,0.0006 ,0.0003 ,-0.0001,-0.0005,-0.0009,-0.0013,-0.0017,-0.0022,-0.0026,-0.0029,-0.0032,-0.0034,-0.0034,-0.0034,-0.0033,-0.0031,-0.0029,-0.0026,-0.0022,-0.0018,-0.0014,-0.0009,-0.0004,0.0001 ,0.0006 ,0.0011 ,0.0015 ,0.0020,
                                                                0.0447,0.0386,0.0330,0.0281,0.0238,0.0201,0.0169,0.0141,0.0118 ,0.0099 ,0.0083 ,0.0069 ,0.0058 ,0.0049 ,0.0041 ,0.0035 ,0.0029 ,0.0025 ,0.0021 ,0.0018 ,0.0015 ,0.0012 ,0.0009 ,0.0005 ,0.0002 ,-0.0002,-0.0006,-0.0010,-0.0014,-0.0018,-0.0022,-0.0026,-0.0028,-0.0030,-0.0031,-0.0031,-0.0030,-0.0028,-0.0026,-0.0023,-0.0020,-0.0016,-0.0012,-0.0008,-0.0003,0.0001 ,0.0006 ,0.0011 ,0.0016 ,0.0020,
                                                                0.0464,0.0403,0.0346,0.0296,0.0251,0.0212,0.0179,0.0150,0.0126 ,0.0106 ,0.0088 ,0.0074 ,0.0062 ,0.0052 ,0.0044 ,0.0038 ,0.0032 ,0.0028 ,0.0024 ,0.0020 ,0.0017 ,0.0014 ,0.0011 ,0.0008 ,0.0005 ,0.0001 ,-0.0003,-0.0007,-0.0011,-0.0015,-0.0019,-0.0022,-0.0025,-0.0026,-0.0027,-0.0027,-0.0027,-0.0025,-0.0023,-0.0021,-0.0018,-0.0014,-0.0010,-0.0007,-0.0002,0.0002 ,0.0006 ,0.0011 ,0.0015 ,0.0020,
                                                                0.0482,0.0419,0.0362,0.0310,0.0265,0.0224,0.0189,0.0160,0.0134 ,0.0112 ,0.0094 ,0.0079 ,0.0067 ,0.0056 ,0.0048 ,0.0041 ,0.0035 ,0.0030 ,0.0026 ,0.0023 ,0.0020 ,0.0017 ,0.0014 ,0.0011 ,0.0008 ,0.0004 ,0.0000 ,-0.0004,-0.0008,-0.0012,-0.0015,-0.0019,-0.0021,-0.0023,-0.0024,-0.0024,-0.0024,-0.0022,-0.0021,-0.0018,-0.0016,-0.0012,-0.0009,-0.0005,-0.0002,0.0002 ,0.0007 ,0.0011 ,0.0015 ,0.0020,
                                                                0.0500,0.0437,0.0378,0.0325,0.0278,0.0236,0.0200,0.0169,0.0142 ,0.0119 ,0.0100 ,0.0084 ,0.0071 ,0.0060 ,0.0051 ,0.0044 ,0.0038 ,0.0033 ,0.0029 ,0.0026 ,0.0023 ,0.0020 ,0.0017 ,0.0014 ,0.0011 ,0.0007 ,0.0004 ,-0.0000,-0.0004,-0.0008,-0.0012,-0.0015,-0.0018,-0.0020,-0.0021,-0.0021,-0.0021,-0.0020,-0.0018,-0.0016,-0.0014,-0.0011,-0.0008,-0.0004,-0.0001,0.0003 ,0.0007 ,0.0011 ,0.0015 ,0.0019,
                                                                0.0518,0.0454,0.0395,0.0340,0.0292,0.0249,0.0211,0.0178,0.0150 ,0.0126 ,0.0106 ,0.0090 ,0.0076 ,0.0064 ,0.0055 ,0.0048 ,0.0041 ,0.0036 ,0.0032 ,0.0029 ,0.0026 ,0.0023 ,0.0020 ,0.0017 ,0.0014 ,0.0011 ,0.0007 ,0.0003 ,-0.0001,-0.0005,-0.0009,-0.0012,-0.0014,-0.0016,-0.0018,-0.0018,-0.0018,-0.0017,-0.0016,-0.0014,-0.0012,-0.0009,-0.0007,-0.0004,-0.0000,0.0003 ,0.0007 ,0.0011 ,0.0015 ,0.0019,
                                                                0.0537,0.0472,0.0411,0.0355,0.0305,0.0261,0.0222,0.0188,0.0158 ,0.0134 ,0.0113 ,0.0095 ,0.0081 ,0.0069 ,0.0059 ,0.0051 ,0.0045 ,0.0040 ,0.0036 ,0.0032 ,0.0029 ,0.0026 ,0.0023 ,0.0021 ,0.0017 ,0.0014 ,0.0010 ,0.0007 ,0.0003 ,-0.0001,-0.0005,-0.0008,-0.0011,-0.0013,-0.0015,-0.0015,-0.0015,-0.0015,-0.0014,-0.0012,-0.0010,-0.0008,-0.0006,-0.0003,0.0000 ,0.0004 ,0.0007 ,0.0011 ,0.0015 ,0.0019,
                                                                0.0555,0.0489,0.0428,0.0371,0.0319,0.0273,0.0233,0.0197,0.0167 ,0.0141 ,0.0119 ,0.0101 ,0.0086 ,0.0073 ,0.0063 ,0.0055 ,0.0048 ,0.0043 ,0.0039 ,0.0035 ,0.0032 ,0.0030 ,0.0027 ,0.0024 ,0.0021 ,0.0018 ,0.0014 ,0.0010 ,0.0006 ,0.0002 ,-0.0002,-0.0005,-0.0008,-0.0010,-0.0012,-0.0013,-0.0013,-0.0012,-0.0012,-0.0010,-0.0009,-0.0007,-0.0005,-0.0002,0.0001 ,0.0004 ,0.0007 ,0.0011 ,0.0014 ,0.0018,
                                                                0.0574,0.0507,0.0444,0.0386,0.0333,0.0286,0.0244,0.0207,0.0175 ,0.0148 ,0.0125 ,0.0106 ,0.0091 ,0.0078 ,0.0067 ,0.0059 ,0.0052 ,0.0047 ,0.0043 ,0.0039 ,0.0036 ,0.0033 ,0.0031 ,0.0028 ,0.0025 ,0.0021 ,0.0018 ,0.0014 ,0.0010 ,0.0005 ,0.0002 ,-0.0002,-0.0005,-0.0007,-0.0009,-0.0010,-0.0011,-0.0010,-0.0010,-0.0009,-0.0007,-0.0006,-0.0004,-0.0001,0.0001 ,0.0004 ,0.0007 ,0.0010 ,0.0014 ,0.0018,
                                                                0.0592,0.0525,0.0461,0.0402,0.0347,0.0298,0.0255,0.0217,0.0184 ,0.0156 ,0.0132 ,0.0112 ,0.0096 ,0.0082 ,0.0072 ,0.0063 ,0.0056 ,0.0051 ,0.0046 ,0.0043 ,0.0040 ,0.0037 ,0.0034 ,0.0032 ,0.0029 ,0.0025 ,0.0021 ,0.0017 ,0.0013 ,0.0009 ,0.0005 ,0.0001 ,-0.0002,-0.0005,-0.0006,-0.0008,-0.0008,-0.0009,-0.0008,-0.0007,-0.0006,-0.0005,-0.0003,-0.0001,0.0002 ,0.0004 ,0.0007 ,0.0010 ,0.0014 ,0.0017,
                                                                0.0611,0.0543,0.0478,0.0417,0.0362,0.0311,0.0266,0.0227,0.0193 ,0.0164 ,0.0139 ,0.0118 ,0.0101 ,0.0087 ,0.0076 ,0.0067 ,0.0060 ,0.0055 ,0.0050 ,0.0047 ,0.0044 ,0.0041 ,0.0038 ,0.0036 ,0.0032 ,0.0029 ,0.0025 ,0.0021 ,0.0016 ,0.0012 ,0.0008 ,0.0004 ,0.0001 ,-0.0002,-0.0004,-0.0005,-0.0006,-0.0007,-0.0007,-0.0006,-0.0005,-0.0004,-0.0002,-0.0000,0.0002 ,0.0004 ,0.0007 ,0.0010 ,0.0013 ,0.0017,
                                                                0.0630,0.0561,0.0495,0.0433,0.0376,0.0324,0.0278,0.0237,0.0202 ,0.0171 ,0.0146 ,0.0124 ,0.0107 ,0.0093 ,0.0081 ,0.0072 ,0.0065 ,0.0059 ,0.0055 ,0.0051 ,0.0048 ,0.0045 ,0.0043 ,0.0040 ,0.0036 ,0.0033 ,0.0029 ,0.0024 ,0.0020 ,0.0015 ,0.0011 ,0.0007 ,0.0003 ,0.0001 ,-0.0002,-0.0003,-0.0005,-0.0005,-0.0005,-0.0005,-0.0004,-0.0003,-0.0002,0.0000 ,0.0002 ,0.0004 ,0.0007 ,0.0010 ,0.0013 ,0.0016,
                                                                0.0649,0.0579,0.0512,0.0449,0.0391,0.0337,0.0290,0.0247,0.0211 ,0.0180 ,0.0153 ,0.0131 ,0.0113 ,0.0098 ,0.0086 ,0.0077 ,0.0069 ,0.0063 ,0.0059 ,0.0055 ,0.0052 ,0.0050 ,0.0047 ,0.0044 ,0.0040 ,0.0037 ,0.0032 ,0.0028 ,0.0023 ,0.0018 ,0.0014 ,0.0010 ,0.0006 ,0.0003 ,0.0000 ,-0.0002,-0.0003,-0.0004,-0.0004,-0.0004,-0.0003,-0.0002,-0.0001,0.0000 ,0.0002 ,0.0004 ,0.0007 ,0.0009 ,0.0012 ,0.0015,
                                                                0.0668,0.0597,0.0530,0.0465,0.0405,0.0351,0.0301,0.0258,0.0220 ,0.0188 ,0.0160 ,0.0138 ,0.0119 ,0.0103 ,0.0091 ,0.0081 ,0.0074 ,0.0068 ,0.0064 ,0.0060 ,0.0057 ,0.0054 ,0.0051 ,0.0048 ,0.0045 ,0.0041 ,0.0036 ,0.0031 ,0.0026 ,0.0021 ,0.0016 ,0.0012 ,0.0008 ,0.0005 ,0.0002 ,0.0000 ,-0.0001,-0.0002,-0.0003,-0.0003,-0.0002,-0.0002,-0.0001,0.0001 ,0.0002 ,0.0004 ,0.0007 ,0.0009 ,0.0012 ,0.0015,
                                                                0.0686,0.0616,0.0547,0.0481,0.0420,0.0364,0.0313,0.0269,0.0230 ,0.0196 ,0.0168 ,0.0144 ,0.0125 ,0.0109 ,0.0097 ,0.0087 ,0.0079 ,0.0073 ,0.0068 ,0.0065 ,0.0062 ,0.0059 ,0.0056 ,0.0052 ,0.0049 ,0.0044 ,0.0040 ,0.0035 ,0.0029 ,0.0024 ,0.0019 ,0.0014 ,0.0010 ,0.0007 ,0.0004 ,0.0002 ,0.0000 ,-0.0001,-0.0002,-0.0002,-0.0002,-0.0001,-0.0000,0.0001 ,0.0003 ,0.0004 ,0.0006 ,0.0009 ,0.0011 ,0.0014,
                                                                0.0705,0.0634,0.0564,0.0497,0.0435,0.0378,0.0326,0.0280,0.0240 ,0.0205 ,0.0176 ,0.0152 ,0.0132 ,0.0115 ,0.0102 ,0.0092 ,0.0084 ,0.0078 ,0.0073 ,0.0070 ,0.0067 ,0.0063 ,0.0060 ,0.0057 ,0.0053 ,0.0048 ,0.0043 ,0.0038 ,0.0032 ,0.0027 ,0.0021 ,0.0017 ,0.0012 ,0.0009 ,0.0006 ,0.0003 ,0.0001 ,0.0000 ,-0.0001,-0.0001,-0.0001,-0.0000,0.0000 ,0.0001 ,0.0003 ,0.0004 ,0.0006 ,0.0008 ,0.0011 ,0.0014,
                                                                0.0724,0.0652,0.0581,0.0514,0.0450,0.0391,0.0338,0.0291,0.0250 ,0.0214 ,0.0184 ,0.0159 ,0.0138 ,0.0122 ,0.0108 ,0.0098 ,0.0090 ,0.0084 ,0.0079 ,0.0075 ,0.0072 ,0.0068 ,0.0065 ,0.0061 ,0.0057 ,0.0052 ,0.0047 ,0.0041 ,0.0035 ,0.0029 ,0.0024 ,0.0019 ,0.0014 ,0.0010 ,0.0007 ,0.0004 ,0.0003 ,0.0001 ,0.0000 ,-0.0000,-0.0000,-0.0000,0.0001 ,0.0001 ,0.0003 ,0.0004 ,0.0006 ,0.0008 ,0.0010 ,0.0013,
                                                                0.0742,0.0670,0.0599,0.0530,0.0465,0.0405,0.0351,0.0302,0.0260 ,0.0224 ,0.0193 ,0.0167 ,0.0145 ,0.0128 ,0.0115 ,0.0104 ,0.0096 ,0.0089 ,0.0084 ,0.0080 ,0.0077 ,0.0073 ,0.0070 ,0.0066 ,0.0061 ,0.0056 ,0.0050 ,0.0044 ,0.0038 ,0.0032 ,0.0026 ,0.0020 ,0.0016 ,0.0012 ,0.0008 ,0.0006 ,0.0004 ,0.0002 ,0.0001 ,0.0000 ,0.0000 ,0.0000 ,0.0001 ,0.0002 ,0.0003 ,0.0004 ,0.0006 ,0.0008 ,0.0010 ,0.0012,
                                                                0.0761,0.0688,0.0616,0.0547,0.0481,0.0419,0.0364,0.0314,0.0271 ,0.0233 ,0.0201 ,0.0175 ,0.0153 ,0.0135 ,0.0121 ,0.0110 ,0.0102 ,0.0095 ,0.0090 ,0.0086 ,0.0082 ,0.0078 ,0.0074 ,0.0070 ,0.0065 ,0.0059 ,0.0053 ,0.0047 ,0.0040 ,0.0034 ,0.0028 ,0.0022 ,0.0017 ,0.0013 ,0.0010 ,0.0007 ,0.0004 ,0.0003 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0001 ,0.0002 ,0.0003 ,0.0004 ,0.0005 ,0.0007 ,0.0009 ,0.0012,
                                                                0.0779,0.0706,0.0634,0.0563,0.0496,0.0434,0.0377,0.0326,0.0282 ,0.0243 ,0.0210 ,0.0183 ,0.0161 ,0.0142 ,0.0128 ,0.0117 ,0.0108 ,0.0101 ,0.0096 ,0.0091 ,0.0087 ,0.0083 ,0.0079 ,0.0074 ,0.0069 ,0.0063 ,0.0056 ,0.0049 ,0.0042 ,0.0036 ,0.0029 ,0.0024 ,0.0019 ,0.0014 ,0.0011 ,0.0008 ,0.0005 ,0.0004 ,0.0002 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0002 ,0.0003 ,0.0004 ,0.0005 ,0.0007 ,0.0009 ,0.0011,
                                                                0.0798,0.0724,0.0651,0.0580,0.0512,0.0448,0.0391,0.0339,0.0293 ,0.0253 ,0.0220 ,0.0192 ,0.0169 ,0.0150 ,0.0135 ,0.0123 ,0.0114 ,0.0107 ,0.0102 ,0.0097 ,0.0093 ,0.0088 ,0.0084 ,0.0078 ,0.0072 ,0.0066 ,0.0059 ,0.0052 ,0.0044 ,0.0038 ,0.0031 ,0.0025 ,0.0020 ,0.0015 ,0.0012 ,0.0008 ,0.0006 ,0.0004 ,0.0003 ,0.0002 ,0.0001 ,0.0001 ,0.0001 ,0.0002 ,0.0003 ,0.0004 ,0.0005 ,0.0007 ,0.0009 ,0.0011,
                                                                0.0816,0.0742,0.0669,0.0597,0.0528,0.0463,0.0404,0.0351,0.0304 ,0.0264 ,0.0229 ,0.0201 ,0.0177 ,0.0158 ,0.0142 ,0.0131 ,0.0121 ,0.0114 ,0.0108 ,0.0103 ,0.0098 ,0.0093 ,0.0088 ,0.0082 ,0.0076 ,0.0069 ,0.0062 ,0.0054 ,0.0046 ,0.0039 ,0.0032 ,0.0026 ,0.0021 ,0.0016 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0004 ,0.0005 ,0.0006 ,0.0008 ,0.0010,
                                                                0.0834,0.0760,0.0686,0.0614,0.0544,0.0478,0.0418,0.0364,0.0316 ,0.0275 ,0.0239 ,0.0210 ,0.0186 ,0.0166 ,0.0150 ,0.0138 ,0.0128 ,0.0121 ,0.0114 ,0.0109 ,0.0104 ,0.0099 ,0.0093 ,0.0087 ,0.0080 ,0.0072 ,0.0064 ,0.0056 ,0.0048 ,0.0041 ,0.0034 ,0.0027 ,0.0022 ,0.0017 ,0.0013 ,0.0010 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0003 ,0.0005 ,0.0006 ,0.0008 ,0.0010,
                                                                0.0852,0.0778,0.0704,0.0631,0.0560,0.0494,0.0433,0.0377,0.0328 ,0.0286 ,0.0250 ,0.0220 ,0.0195 ,0.0174 ,0.0158 ,0.0145 ,0.0135 ,0.0127 ,0.0121 ,0.0115 ,0.0109 ,0.0104 ,0.0097 ,0.0090 ,0.0083 ,0.0075 ,0.0066 ,0.0058 ,0.0050 ,0.0042 ,0.0035 ,0.0028 ,0.0023 ,0.0018 ,0.0014 ,0.0010 ,0.0008 ,0.0006 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0004 ,0.0006 ,0.0007 ,0.0009,
                                                                0.0870,0.0796,0.0721,0.0648,0.0576,0.0509,0.0447,0.0391,0.0341 ,0.0298 ,0.0261 ,0.0230 ,0.0204 ,0.0183 ,0.0167 ,0.0153 ,0.0143 ,0.0135 ,0.0127 ,0.0121 ,0.0115 ,0.0109 ,0.0102 ,0.0094 ,0.0086 ,0.0077 ,0.0069 ,0.0060 ,0.0051 ,0.0043 ,0.0036 ,0.0029 ,0.0023 ,0.0018 ,0.0014 ,0.0011 ,0.0008 ,0.0006 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0004 ,0.0005 ,0.0007 ,0.0009,
                                                                0.0887,0.0814,0.0739,0.0665,0.0593,0.0525,0.0462,0.0405,0.0354 ,0.0310 ,0.0272 ,0.0240 ,0.0214 ,0.0192 ,0.0175 ,0.0162 ,0.0151 ,0.0142 ,0.0134 ,0.0127 ,0.0120 ,0.0113 ,0.0106 ,0.0098 ,0.0089 ,0.0080 ,0.0071 ,0.0061 ,0.0052 ,0.0044 ,0.0037 ,0.0030 ,0.0024 ,0.0019 ,0.0015 ,0.0011 ,0.0008 ,0.0006 ,0.0005 ,0.0003 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0004 ,0.0005 ,0.0006 ,0.0008,
                                                                0.0905,0.0832,0.0757,0.0682,0.0610,0.0541,0.0477,0.0419,0.0367 ,0.0322 ,0.0283 ,0.0251 ,0.0224 ,0.0202 ,0.0184 ,0.0170 ,0.0158 ,0.0149 ,0.0141 ,0.0133 ,0.0126 ,0.0118 ,0.0110 ,0.0101 ,0.0092 ,0.0082 ,0.0072 ,0.0063 ,0.0054 ,0.0045 ,0.0037 ,0.0030 ,0.0024 ,0.0019 ,0.0015 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0004 ,0.0005 ,0.0006 ,0.0008,
                                                                0.0922,0.0849,0.0774,0.0700,0.0627,0.0557,0.0493,0.0434,0.0381 ,0.0335 ,0.0295 ,0.0262 ,0.0234 ,0.0212 ,0.0193 ,0.0179 ,0.0167 ,0.0156 ,0.0148 ,0.0139 ,0.0131 ,0.0123 ,0.0114 ,0.0104 ,0.0094 ,0.0084 ,0.0074 ,0.0064 ,0.0054 ,0.0046 ,0.0038 ,0.0031 ,0.0025 ,0.0020 ,0.0015 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0003 ,0.0004 ,0.0006 ,0.0007,
                                                                0.0940,0.0867,0.0792,0.0717,0.0644,0.0574,0.0509,0.0449,0.0395 ,0.0348 ,0.0308 ,0.0274 ,0.0245 ,0.0222 ,0.0203 ,0.0187 ,0.0175 ,0.0164 ,0.0154 ,0.0146 ,0.0137 ,0.0127 ,0.0118 ,0.0107 ,0.0097 ,0.0086 ,0.0075 ,0.0065 ,0.0055 ,0.0046 ,0.0038 ,0.0031 ,0.0025 ,0.0020 ,0.0016 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0004 ,0.0005 ,0.0007,
                                                                0.0957,0.0884,0.0810,0.0735,0.0661,0.0591,0.0525,0.0464,0.0410 ,0.0362 ,0.0321 ,0.0285 ,0.0256 ,0.0232 ,0.0213 ,0.0196 ,0.0183 ,0.0172 ,0.0161 ,0.0152 ,0.0142 ,0.0132 ,0.0121 ,0.0110 ,0.0099 ,0.0088 ,0.0077 ,0.0066 ,0.0056 ,0.0047 ,0.0039 ,0.0031 ,0.0025 ,0.0020 ,0.0016 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0004 ,0.0005 ,0.0006,
                                                                0.0974,0.0902,0.0828,0.0753,0.0679,0.0608,0.0541,0.0480,0.0425 ,0.0376 ,0.0334 ,0.0298 ,0.0268 ,0.0243 ,0.0223 ,0.0206 ,0.0192 ,0.0179 ,0.0168 ,0.0157 ,0.0147 ,0.0136 ,0.0125 ,0.0113 ,0.0101 ,0.0089 ,0.0078 ,0.0067 ,0.0056 ,0.0047 ,0.0039 ,0.0032 ,0.0025 ,0.0020 ,0.0016 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0003 ,0.0004 ,0.0006,
                                                                0.0990,0.0919,0.0845,0.0771,0.0697,0.0626,0.0558,0.0496,0.0440 ,0.0390 ,0.0347 ,0.0311 ,0.0280 ,0.0254 ,0.0233 ,0.0215 ,0.0200 ,0.0187 ,0.0175 ,0.0163 ,0.0152 ,0.0140 ,0.0128 ,0.0115 ,0.0103 ,0.0090 ,0.0078 ,0.0067 ,0.0057 ,0.0047 ,0.0039 ,0.0032 ,0.0025 ,0.0020 ,0.0016 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0004 ,0.0005,
                                                                0.1007,0.0937,0.0863,0.0789,0.0715,0.0643,0.0575,0.0513,0.0456 ,0.0405 ,0.0361 ,0.0324 ,0.0292 ,0.0266 ,0.0243 ,0.0225 ,0.0209 ,0.0195 ,0.0182 ,0.0169 ,0.0157 ,0.0144 ,0.0131 ,0.0118 ,0.0105 ,0.0092 ,0.0079 ,0.0068 ,0.0057 ,0.0047 ,0.0039 ,0.0032 ,0.0025 ,0.0020 ,0.0016 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0004 ,0.0005,
                                                                0.1024,0.0954,0.0881,0.0807,0.0733,0.0661,0.0593,0.0530,0.0472 ,0.0421 ,0.0376 ,0.0337 ,0.0305 ,0.0277 ,0.0254 ,0.0234 ,0.0218 ,0.0202 ,0.0188 ,0.0175 ,0.0161 ,0.0147 ,0.0134 ,0.0120 ,0.0106 ,0.0093 ,0.0080 ,0.0068 ,0.0057 ,0.0047 ,0.0039 ,0.0032 ,0.0025 ,0.0020 ,0.0016 ,0.0012 ,0.0010 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0002 ,0.0003 ,0.0004 ,0.0005};
        
        return arrayCurvatureMin[y * arrayWidth + x];
    }
    case 10:
    {
        static const std::vector<double> arrayCurvatureMin   = {0.0270,0.0206,0.0152,0.0107,0.0072,0.0045,0.0025,0.0012,0.0004,-0.0000,-0.0002,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0004,-0.0006,-0.0007,-0.0009,-0.0011,-0.0013,-0.0015,-0.0017,-0.0019,-0.0021,-0.0022,-0.0024,-0.0026,-0.0027,-0.0029,-0.0032,-0.0034,-0.0037,-0.0041,-0.0044,-0.0048,-0.0051,-0.0054,-0.0056,-0.0056,-0.0055,-0.0053,-0.0050,-0.0046,-0.0042,-0.0037,-0.0032,-0.0028,-0.0025,
                                                                0.0285,0.0220,0.0165,0.0119,0.0082,0.0054,0.0034,0.0019,0.0010,0.0005 ,0.0003 ,0.0001 ,0.0001 ,0.0001 ,0.0000 ,-0.0000,-0.0002,-0.0003,-0.0005,-0.0007,-0.0009,-0.0011,-0.0014,-0.0016,-0.0018,-0.0019,-0.0021,-0.0023,-0.0024,-0.0026,-0.0028,-0.0030,-0.0033,-0.0036,-0.0039,-0.0043,-0.0046,-0.0050,-0.0052,-0.0054,-0.0054,-0.0053,-0.0051,-0.0048,-0.0044,-0.0039,-0.0035,-0.0030,-0.0026,-0.0023,
                                                                0.0300,0.0234,0.0178,0.0130,0.0093,0.0063,0.0042,0.0027,0.0017,0.0011 ,0.0008 ,0.0006 ,0.0005 ,0.0004 ,0.0004 ,0.0002 ,0.0001 ,-0.0001,-0.0003,-0.0005,-0.0007,-0.0010,-0.0012,-0.0014,-0.0016,-0.0018,-0.0020,-0.0021,-0.0023,-0.0024,-0.0026,-0.0029,-0.0031,-0.0034,-0.0038,-0.0041,-0.0045,-0.0048,-0.0050,-0.0052,-0.0052,-0.0051,-0.0048,-0.0045,-0.0041,-0.0037,-0.0033,-0.0028,-0.0025,-0.0022,
                                                                0.0315,0.0248,0.0191,0.0142,0.0103,0.0073,0.0050,0.0034,0.0024,0.0017 ,0.0013 ,0.0010 ,0.0009 ,0.0008 ,0.0007 ,0.0005 ,0.0004 ,0.0002 ,-0.0000,-0.0003,-0.0005,-0.0008,-0.0010,-0.0012,-0.0014,-0.0016,-0.0018,-0.0020,-0.0021,-0.0023,-0.0025,-0.0027,-0.0030,-0.0033,-0.0036,-0.0039,-0.0043,-0.0046,-0.0048,-0.0049,-0.0049,-0.0048,-0.0046,-0.0043,-0.0039,-0.0035,-0.0030,-0.0026,-0.0023,-0.0020,
                                                                0.0330,0.0262,0.0204,0.0154,0.0114,0.0083,0.0059,0.0042,0.0030,0.0023 ,0.0018 ,0.0015 ,0.0013 ,0.0012 ,0.0010 ,0.0009 ,0.0007 ,0.0004 ,0.0002 ,-0.0001,-0.0003,-0.0006,-0.0008,-0.0011,-0.0013,-0.0015,-0.0016,-0.0018,-0.0020,-0.0021,-0.0023,-0.0025,-0.0028,-0.0031,-0.0034,-0.0038,-0.0041,-0.0044,-0.0046,-0.0047,-0.0047,-0.0046,-0.0044,-0.0040,-0.0037,-0.0032,-0.0028,-0.0025,-0.0021,-0.0019,
                                                                0.0345,0.0277,0.0217,0.0167,0.0125,0.0093,0.0068,0.0050,0.0037,0.0029 ,0.0024 ,0.0020 ,0.0018 ,0.0016 ,0.0014 ,0.0012 ,0.0009 ,0.0007 ,0.0004 ,0.0002 ,-0.0001,-0.0004,-0.0007,-0.0009,-0.0011,-0.0013,-0.0015,-0.0016,-0.0018,-0.0020,-0.0022,-0.0024,-0.0026,-0.0029,-0.0032,-0.0036,-0.0039,-0.0042,-0.0044,-0.0045,-0.0045,-0.0043,-0.0041,-0.0038,-0.0034,-0.0030,-0.0026,-0.0023,-0.0020,-0.0017,
                                                                0.0361,0.0291,0.0231,0.0179,0.0137,0.0103,0.0077,0.0058,0.0045,0.0035 ,0.0029 ,0.0025 ,0.0022 ,0.0020 ,0.0017 ,0.0015 ,0.0012 ,0.0010 ,0.0007 ,0.0004 ,0.0001 ,-0.0002,-0.0005,-0.0007,-0.0009,-0.0011,-0.0013,-0.0015,-0.0016,-0.0018,-0.0020,-0.0022,-0.0025,-0.0028,-0.0031,-0.0034,-0.0037,-0.0040,-0.0042,-0.0043,-0.0042,-0.0041,-0.0039,-0.0036,-0.0032,-0.0028,-0.0024,-0.0021,-0.0018,-0.0016,
                                                                0.0376,0.0306,0.0244,0.0191,0.0148,0.0113,0.0086,0.0066,0.0052,0.0042 ,0.0035 ,0.0030 ,0.0027 ,0.0024 ,0.0021 ,0.0018 ,0.0015 ,0.0012 ,0.0009 ,0.0006 ,0.0003 ,0.0000 ,-0.0003,-0.0005,-0.0007,-0.0009,-0.0011,-0.0013,-0.0015,-0.0016,-0.0018,-0.0020,-0.0023,-0.0026,-0.0029,-0.0032,-0.0035,-0.0038,-0.0040,-0.0040,-0.0040,-0.0039,-0.0036,-0.0033,-0.0030,-0.0026,-0.0022,-0.0019,-0.0017,-0.0015,
                                                                0.0392,0.0321,0.0258,0.0204,0.0160,0.0124,0.0096,0.0075,0.0059,0.0048 ,0.0041 ,0.0035 ,0.0031 ,0.0028 ,0.0025 ,0.0022 ,0.0019 ,0.0015 ,0.0012 ,0.0009 ,0.0005 ,0.0002 ,-0.0001,-0.0003,-0.0006,-0.0008,-0.0009,-0.0011,-0.0013,-0.0014,-0.0016,-0.0018,-0.0021,-0.0024,-0.0027,-0.0030,-0.0033,-0.0036,-0.0037,-0.0038,-0.0038,-0.0036,-0.0034,-0.0031,-0.0028,-0.0024,-0.0021,-0.0018,-0.0015,-0.0014,
                                                                0.0408,0.0335,0.0272,0.0217,0.0171,0.0134,0.0105,0.0083,0.0067,0.0055 ,0.0047 ,0.0041 ,0.0036 ,0.0032 ,0.0029 ,0.0025 ,0.0022 ,0.0018 ,0.0015 ,0.0011 ,0.0008 ,0.0004 ,0.0001 ,-0.0001,-0.0004,-0.0006,-0.0008,-0.0009,-0.0011,-0.0013,-0.0014,-0.0017,-0.0019,-0.0022,-0.0025,-0.0028,-0.0031,-0.0033,-0.0035,-0.0036,-0.0035,-0.0034,-0.0032,-0.0029,-0.0025,-0.0022,-0.0019,-0.0016,-0.0014,-0.0013,
                                                                0.0423,0.0350,0.0286,0.0230,0.0183,0.0145,0.0115,0.0092,0.0075,0.0062 ,0.0053 ,0.0046 ,0.0041 ,0.0037 ,0.0033 ,0.0029 ,0.0025 ,0.0021 ,0.0017 ,0.0014 ,0.0010 ,0.0007 ,0.0004 ,0.0001 ,-0.0002,-0.0004,-0.0006,-0.0007,-0.0009,-0.0011,-0.0013,-0.0015,-0.0017,-0.0020,-0.0023,-0.0026,-0.0029,-0.0031,-0.0033,-0.0033,-0.0033,-0.0032,-0.0029,-0.0027,-0.0023,-0.0020,-0.0017,-0.0014,-0.0013,-0.0011,
                                                                0.0439,0.0366,0.0300,0.0243,0.0195,0.0156,0.0125,0.0101,0.0083,0.0069 ,0.0059 ,0.0052 ,0.0046 ,0.0041 ,0.0037 ,0.0033 ,0.0029 ,0.0024 ,0.0020 ,0.0016 ,0.0013 ,0.0009 ,0.0006 ,0.0003 ,0.0000 ,-0.0002,-0.0004,-0.0005,-0.0007,-0.0009,-0.0011,-0.0013,-0.0015,-0.0018,-0.0021,-0.0024,-0.0027,-0.0029,-0.0031,-0.0031,-0.0031,-0.0029,-0.0027,-0.0024,-0.0021,-0.0018,-0.0015,-0.0013,-0.0011,-0.0010,
                                                                0.0455,0.0381,0.0314,0.0256,0.0207,0.0167,0.0135,0.0110,0.0091,0.0076 ,0.0066 ,0.0058 ,0.0051 ,0.0046 ,0.0041 ,0.0036 ,0.0032 ,0.0028 ,0.0023 ,0.0019 ,0.0015 ,0.0011 ,0.0008 ,0.0005 ,0.0003 ,0.0000 ,-0.0002,-0.0003,-0.0005,-0.0007,-0.0009,-0.0011,-0.0013,-0.0016,-0.0019,-0.0022,-0.0025,-0.0027,-0.0028,-0.0029,-0.0028,-0.0027,-0.0025,-0.0022,-0.0019,-0.0016,-0.0014,-0.0012,-0.0010,-0.0009,
                                                                0.0471,0.0396,0.0329,0.0270,0.0220,0.0178,0.0145,0.0119,0.0099,0.0084 ,0.0072 ,0.0063 ,0.0057 ,0.0051 ,0.0045 ,0.0040 ,0.0036 ,0.0031 ,0.0026 ,0.0022 ,0.0018 ,0.0014 ,0.0010 ,0.0007 ,0.0005 ,0.0002 ,0.0000 ,-0.0001,-0.0003,-0.0005,-0.0006,-0.0009,-0.0011,-0.0014,-0.0017,-0.0020,-0.0023,-0.0025,-0.0026,-0.0027,-0.0026,-0.0025,-0.0023,-0.0020,-0.0017,-0.0015,-0.0012,-0.0010,-0.0009,-0.0008,
                                                                0.0487,0.0411,0.0343,0.0283,0.0232,0.0190,0.0155,0.0128,0.0107,0.0091 ,0.0079 ,0.0070 ,0.0062 ,0.0055 ,0.0050 ,0.0044 ,0.0039 ,0.0034 ,0.0029 ,0.0025 ,0.0020 ,0.0016 ,0.0013 ,0.0010 ,0.0007 ,0.0005 ,0.0003 ,0.0001 ,-0.0001,-0.0002,-0.0004,-0.0006,-0.0009,-0.0012,-0.0015,-0.0018,-0.0020,-0.0022,-0.0024,-0.0024,-0.0024,-0.0023,-0.0021,-0.0018,-0.0015,-0.0013,-0.0011,-0.0009,-0.0008,-0.0008,
                                                                0.0503,0.0427,0.0358,0.0297,0.0245,0.0201,0.0166,0.0138,0.0116,0.0099 ,0.0086 ,0.0076 ,0.0067 ,0.0060 ,0.0054 ,0.0049 ,0.0043 ,0.0038 ,0.0033 ,0.0028 ,0.0023 ,0.0019 ,0.0015 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0003 ,0.0001 ,-0.0000,-0.0002,-0.0004,-0.0007,-0.0009,-0.0012,-0.0015,-0.0018,-0.0020,-0.0022,-0.0022,-0.0022,-0.0020,-0.0018,-0.0016,-0.0014,-0.0011,-0.0009,-0.0008,-0.0007,-0.0007,
                                                                0.0520,0.0442,0.0373,0.0311,0.0258,0.0213,0.0177,0.0147,0.0124,0.0107 ,0.0093 ,0.0082 ,0.0073 ,0.0066 ,0.0059 ,0.0053 ,0.0047 ,0.0041 ,0.0036 ,0.0031 ,0.0026 ,0.0022 ,0.0018 ,0.0015 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0002 ,0.0000 ,-0.0002,-0.0004,-0.0007,-0.0010,-0.0013,-0.0016,-0.0018,-0.0019,-0.0020,-0.0019,-0.0018,-0.0016,-0.0014,-0.0012,-0.0010,-0.0008,-0.0006,-0.0006,-0.0006,
                                                                0.0536,0.0458,0.0388,0.0325,0.0271,0.0225,0.0187,0.0157,0.0133,0.0115 ,0.0100 ,0.0088 ,0.0079 ,0.0071 ,0.0064 ,0.0057 ,0.0051 ,0.0045 ,0.0039 ,0.0034 ,0.0029 ,0.0025 ,0.0021 ,0.0017 ,0.0014 ,0.0012 ,0.0010 ,0.0008 ,0.0006 ,0.0004 ,0.0002 ,0.0000 ,-0.0002,-0.0005,-0.0008,-0.0011,-0.0014,-0.0016,-0.0017,-0.0018,-0.0017,-0.0016,-0.0014,-0.0012,-0.0010,-0.0008,-0.0006,-0.0005,-0.0005,-0.0005,
                                                                0.0552,0.0474,0.0403,0.0339,0.0284,0.0237,0.0199,0.0167,0.0142,0.0123 ,0.0107 ,0.0095 ,0.0085 ,0.0076 ,0.0069 ,0.0062 ,0.0055 ,0.0049 ,0.0043 ,0.0037 ,0.0032 ,0.0027 ,0.0023 ,0.0020 ,0.0017 ,0.0014 ,0.0012 ,0.0010 ,0.0009 ,0.0007 ,0.0005 ,0.0003 ,0.0000 ,-0.0003,-0.0006,-0.0009,-0.0011,-0.0013,-0.0015,-0.0015,-0.0015,-0.0014,-0.0012,-0.0010,-0.0008,-0.0006,-0.0005,-0.0004,-0.0004,-0.0004,
                                                                0.0569,0.0490,0.0418,0.0353,0.0297,0.0249,0.0210,0.0177,0.0152,0.0131 ,0.0115 ,0.0102 ,0.0091 ,0.0082 ,0.0074 ,0.0066 ,0.0059 ,0.0052 ,0.0046 ,0.0040 ,0.0035 ,0.0030 ,0.0026 ,0.0023 ,0.0019 ,0.0017 ,0.0015 ,0.0013 ,0.0011 ,0.0009 ,0.0007 ,0.0005 ,0.0003 ,-0.0000,-0.0003,-0.0006,-0.0009,-0.0011,-0.0013,-0.0013,-0.0013,-0.0012,-0.0010,-0.0009,-0.0007,-0.0005,-0.0004,-0.0003,-0.0003,-0.0004,
                                                                0.0585,0.0506,0.0433,0.0368,0.0311,0.0262,0.0221,0.0188,0.0161,0.0139 ,0.0122 ,0.0109 ,0.0097 ,0.0087 ,0.0079 ,0.0071 ,0.0063 ,0.0056 ,0.0050 ,0.0044 ,0.0038 ,0.0033 ,0.0029 ,0.0025 ,0.0022 ,0.0019 ,0.0017 ,0.0015 ,0.0014 ,0.0012 ,0.0010 ,0.0008 ,0.0005 ,0.0002 ,-0.0001,-0.0004,-0.0007,-0.0009,-0.0010,-0.0011,-0.0011,-0.0010,-0.0008,-0.0007,-0.0005,-0.0004,-0.0002,-0.0002,-0.0002,-0.0003,
                                                                0.0602,0.0522,0.0448,0.0382,0.0324,0.0274,0.0233,0.0198,0.0170,0.0148 ,0.0130 ,0.0116 ,0.0104 ,0.0093 ,0.0084 ,0.0076 ,0.0068 ,0.0060 ,0.0054 ,0.0047 ,0.0042 ,0.0037 ,0.0032 ,0.0028 ,0.0025 ,0.0022 ,0.0020 ,0.0018 ,0.0016 ,0.0014 ,0.0012 ,0.0010 ,0.0007 ,0.0004 ,0.0001 ,-0.0002,-0.0005,-0.0007,-0.0008,-0.0009,-0.0009,-0.0008,-0.0007,-0.0005,-0.0004,-0.0002,-0.0001,-0.0001,-0.0001,-0.0003,
                                                                0.0618,0.0538,0.0464,0.0397,0.0338,0.0287,0.0244,0.0209,0.0180,0.0157 ,0.0138 ,0.0123 ,0.0110 ,0.0099 ,0.0089 ,0.0080 ,0.0072 ,0.0065 ,0.0057 ,0.0051 ,0.0045 ,0.0040 ,0.0035 ,0.0031 ,0.0028 ,0.0025 ,0.0023 ,0.0021 ,0.0019 ,0.0017 ,0.0015 ,0.0013 ,0.0010 ,0.0007 ,0.0004 ,0.0001 ,-0.0002,-0.0004,-0.0006,-0.0007,-0.0006,-0.0006,-0.0005,-0.0003,-0.0002,-0.0001,-0.0000,-0.0000,-0.0001,-0.0002,
                                                                0.0635,0.0554,0.0479,0.0412,0.0352,0.0300,0.0256,0.0220,0.0190,0.0166 ,0.0146 ,0.0130 ,0.0117 ,0.0105 ,0.0095 ,0.0085 ,0.0077 ,0.0069 ,0.0061 ,0.0055 ,0.0048 ,0.0043 ,0.0038 ,0.0034 ,0.0031 ,0.0028 ,0.0026 ,0.0024 ,0.0022 ,0.0020 ,0.0018 ,0.0015 ,0.0012 ,0.0009 ,0.0006 ,0.0003 ,0.0000 ,-0.0002,-0.0004,-0.0004,-0.0004,-0.0004,-0.0003,-0.0002,-0.0001,0.0000 ,0.0001 ,0.0001 ,0.0000 ,-0.0002,
                                                                0.0651,0.0570,0.0495,0.0426,0.0366,0.0313,0.0268,0.0231,0.0200,0.0175 ,0.0154 ,0.0138 ,0.0123 ,0.0111 ,0.0100 ,0.0091 ,0.0082 ,0.0073 ,0.0066 ,0.0058 ,0.0052 ,0.0046 ,0.0041 ,0.0037 ,0.0034 ,0.0031 ,0.0028 ,0.0026 ,0.0024 ,0.0023 ,0.0020 ,0.0018 ,0.0015 ,0.0012 ,0.0008 ,0.0005 ,0.0002 ,0.0000 ,-0.0002,-0.0002,-0.0002,-0.0002,-0.0001,-0.0000,0.0001 ,0.0002 ,0.0002 ,0.0002 ,0.0001 ,-0.0001,
                                                                0.0668,0.0586,0.0510,0.0441,0.0380,0.0326,0.0280,0.0242,0.0210,0.0184 ,0.0163 ,0.0145 ,0.0130 ,0.0118 ,0.0106 ,0.0096 ,0.0086 ,0.0078 ,0.0070 ,0.0062 ,0.0056 ,0.0050 ,0.0045 ,0.0041 ,0.0037 ,0.0034 ,0.0031 ,0.0029 ,0.0027 ,0.0025 ,0.0023 ,0.0020 ,0.0018 ,0.0014 ,0.0011 ,0.0008 ,0.0005 ,0.0002 ,0.0001 ,-0.0000,-0.0001,-0.0000,0.0001 ,0.0001 ,0.0002 ,0.0003 ,0.0003 ,0.0002 ,0.0001 ,-0.0001,
                                                                0.0684,0.0602,0.0526,0.0456,0.0394,0.0339,0.0293,0.0253,0.0221,0.0194 ,0.0171 ,0.0153 ,0.0137 ,0.0124 ,0.0112 ,0.0101 ,0.0091 ,0.0082 ,0.0074 ,0.0066 ,0.0060 ,0.0053 ,0.0048 ,0.0044 ,0.0040 ,0.0037 ,0.0035 ,0.0032 ,0.0030 ,0.0028 ,0.0026 ,0.0023 ,0.0020 ,0.0017 ,0.0013 ,0.0010 ,0.0007 ,0.0004 ,0.0003 ,0.0002 ,0.0001 ,0.0002 ,0.0002 ,0.0003 ,0.0004 ,0.0004 ,0.0004 ,0.0003 ,0.0002 ,-0.0000,
                                                                0.0701,0.0619,0.0542,0.0472,0.0408,0.0353,0.0305,0.0265,0.0231,0.0203 ,0.0180 ,0.0161 ,0.0145 ,0.0131 ,0.0118 ,0.0107 ,0.0096 ,0.0087 ,0.0078 ,0.0070 ,0.0063 ,0.0057 ,0.0052 ,0.0047 ,0.0043 ,0.0040 ,0.0038 ,0.0035 ,0.0033 ,0.0031 ,0.0029 ,0.0026 ,0.0023 ,0.0019 ,0.0016 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0004 ,0.0003 ,0.0003 ,0.0004 ,0.0004 ,0.0005 ,0.0005 ,0.0005 ,0.0004 ,0.0002 ,-0.0000,
                                                                0.0718,0.0635,0.0558,0.0487,0.0423,0.0366,0.0318,0.0276,0.0242,0.0213 ,0.0189 ,0.0169 ,0.0152 ,0.0137 ,0.0124 ,0.0112 ,0.0102 ,0.0092 ,0.0083 ,0.0075 ,0.0067 ,0.0061 ,0.0055 ,0.0051 ,0.0047 ,0.0044 ,0.0041 ,0.0039 ,0.0036 ,0.0034 ,0.0032 ,0.0029 ,0.0025 ,0.0022 ,0.0018 ,0.0015 ,0.0011 ,0.0009 ,0.0007 ,0.0006 ,0.0005 ,0.0005 ,0.0005 ,0.0006 ,0.0006 ,0.0006 ,0.0006 ,0.0005 ,0.0003 ,0.0000 ,
                                                                0.0734,0.0652,0.0574,0.0502,0.0437,0.0380,0.0331,0.0288,0.0253,0.0223 ,0.0198 ,0.0177 ,0.0160 ,0.0144 ,0.0131 ,0.0118 ,0.0107 ,0.0097 ,0.0087 ,0.0079 ,0.0071 ,0.0065 ,0.0059 ,0.0054 ,0.0050 ,0.0047 ,0.0044 ,0.0042 ,0.0040 ,0.0037 ,0.0035 ,0.0032 ,0.0028 ,0.0024 ,0.0021 ,0.0017 ,0.0014 ,0.0011 ,0.0009 ,0.0008 ,0.0007 ,0.0007 ,0.0007 ,0.0007 ,0.0007 ,0.0007 ,0.0006 ,0.0005 ,0.0003 ,0.0000 ,
                                                                0.0751,0.0668,0.0590,0.0518,0.0452,0.0394,0.0343,0.0300,0.0264,0.0233 ,0.0207 ,0.0186 ,0.0167 ,0.0151 ,0.0137 ,0.0124 ,0.0113 ,0.0102 ,0.0092 ,0.0083 ,0.0076 ,0.0069 ,0.0063 ,0.0058 ,0.0054 ,0.0051 ,0.0048 ,0.0045 ,0.0043 ,0.0040 ,0.0038 ,0.0034 ,0.0031 ,0.0027 ,0.0023 ,0.0019 ,0.0016 ,0.0013 ,0.0011 ,0.0010 ,0.0009 ,0.0008 ,0.0008 ,0.0008 ,0.0008 ,0.0008 ,0.0007 ,0.0006 ,0.0004 ,0.0001 ,
                                                                0.0767,0.0684,0.0606,0.0533,0.0467,0.0408,0.0357,0.0312,0.0275,0.0243 ,0.0217 ,0.0194 ,0.0175 ,0.0158 ,0.0144 ,0.0130 ,0.0118 ,0.0107 ,0.0097 ,0.0088 ,0.0080 ,0.0073 ,0.0067 ,0.0062 ,0.0058 ,0.0054 ,0.0051 ,0.0049 ,0.0046 ,0.0043 ,0.0041 ,0.0037 ,0.0034 ,0.0030 ,0.0025 ,0.0022 ,0.0018 ,0.0015 ,0.0013 ,0.0011 ,0.0010 ,0.0010 ,0.0010 ,0.0010 ,0.0009 ,0.0009 ,0.0008 ,0.0006 ,0.0004 ,0.0001 ,
                                                                0.0784,0.0701,0.0622,0.0549,0.0482,0.0422,0.0370,0.0325,0.0286,0.0254 ,0.0226 ,0.0203 ,0.0183 ,0.0166 ,0.0150 ,0.0137 ,0.0124 ,0.0112 ,0.0102 ,0.0093 ,0.0084 ,0.0077 ,0.0071 ,0.0066 ,0.0061 ,0.0058 ,0.0055 ,0.0052 ,0.0049 ,0.0047 ,0.0044 ,0.0040 ,0.0036 ,0.0032 ,0.0028 ,0.0024 ,0.0020 ,0.0017 ,0.0015 ,0.0013 ,0.0012 ,0.0011 ,0.0011 ,0.0011 ,0.0010 ,0.0010 ,0.0009 ,0.0007 ,0.0004 ,0.0001 ,
                                                                0.0800,0.0717,0.0638,0.0564,0.0497,0.0436,0.0383,0.0337,0.0298,0.0264 ,0.0236 ,0.0212 ,0.0191 ,0.0173 ,0.0157 ,0.0143 ,0.0130 ,0.0118 ,0.0107 ,0.0097 ,0.0089 ,0.0081 ,0.0075 ,0.0070 ,0.0065 ,0.0062 ,0.0058 ,0.0056 ,0.0053 ,0.0050 ,0.0047 ,0.0043 ,0.0039 ,0.0035 ,0.0030 ,0.0026 ,0.0022 ,0.0019 ,0.0017 ,0.0015 ,0.0014 ,0.0013 ,0.0012 ,0.0012 ,0.0011 ,0.0011 ,0.0009 ,0.0007 ,0.0005 ,0.0001 ,
                                                                0.0817,0.0734,0.0654,0.0580,0.0512,0.0451,0.0397,0.0350,0.0309,0.0275 ,0.0246 ,0.0221 ,0.0200 ,0.0181 ,0.0164 ,0.0149 ,0.0136 ,0.0124 ,0.0112 ,0.0102 ,0.0094 ,0.0086 ,0.0079 ,0.0074 ,0.0069 ,0.0065 ,0.0062 ,0.0059 ,0.0056 ,0.0053 ,0.0050 ,0.0046 ,0.0042 ,0.0037 ,0.0033 ,0.0028 ,0.0024 ,0.0021 ,0.0019 ,0.0017 ,0.0015 ,0.0014 ,0.0014 ,0.0013 ,0.0012 ,0.0011 ,0.0010 ,0.0008 ,0.0005 ,0.0001 ,
                                                                0.0834,0.0750,0.0671,0.0596,0.0527,0.0465,0.0410,0.0362,0.0321,0.0286 ,0.0256 ,0.0230 ,0.0208 ,0.0189 ,0.0172 ,0.0156 ,0.0142 ,0.0129 ,0.0118 ,0.0107 ,0.0098 ,0.0090 ,0.0084 ,0.0078 ,0.0073 ,0.0069 ,0.0066 ,0.0063 ,0.0060 ,0.0057 ,0.0053 ,0.0049 ,0.0045 ,0.0040 ,0.0035 ,0.0031 ,0.0027 ,0.0023 ,0.0020 ,0.0018 ,0.0017 ,0.0016 ,0.0015 ,0.0014 ,0.0013 ,0.0012 ,0.0010 ,0.0008 ,0.0005 ,0.0001 ,
                                                                0.0850,0.0767,0.0687,0.0612,0.0542,0.0480,0.0424,0.0375,0.0333,0.0297 ,0.0266 ,0.0240 ,0.0217 ,0.0197 ,0.0179 ,0.0163 ,0.0148 ,0.0135 ,0.0123 ,0.0113 ,0.0103 ,0.0095 ,0.0088 ,0.0082 ,0.0077 ,0.0073 ,0.0070 ,0.0067 ,0.0063 ,0.0060 ,0.0056 ,0.0052 ,0.0047 ,0.0042 ,0.0038 ,0.0033 ,0.0029 ,0.0025 ,0.0022 ,0.0020 ,0.0018 ,0.0017 ,0.0016 ,0.0015 ,0.0014 ,0.0013 ,0.0011 ,0.0009 ,0.0005 ,0.0001 ,
                                                                0.0866,0.0783,0.0703,0.0628,0.0558,0.0494,0.0438,0.0388,0.0345,0.0308 ,0.0277 ,0.0249 ,0.0226 ,0.0205 ,0.0186 ,0.0170 ,0.0155 ,0.0141 ,0.0129 ,0.0118 ,0.0108 ,0.0100 ,0.0093 ,0.0087 ,0.0082 ,0.0078 ,0.0074 ,0.0070 ,0.0067 ,0.0063 ,0.0059 ,0.0055 ,0.0050 ,0.0045 ,0.0040 ,0.0035 ,0.0031 ,0.0027 ,0.0024 ,0.0022 ,0.0020 ,0.0019 ,0.0017 ,0.0016 ,0.0015 ,0.0014 ,0.0012 ,0.0009 ,0.0005 ,0.0001 ,
                                                                0.0883,0.0800,0.0719,0.0644,0.0573,0.0509,0.0452,0.0402,0.0358,0.0320 ,0.0287 ,0.0259 ,0.0235 ,0.0213 ,0.0194 ,0.0177 ,0.0161 ,0.0147 ,0.0135 ,0.0123 ,0.0114 ,0.0105 ,0.0098 ,0.0091 ,0.0086 ,0.0082 ,0.0078 ,0.0074 ,0.0071 ,0.0067 ,0.0063 ,0.0058 ,0.0053 ,0.0048 ,0.0042 ,0.0037 ,0.0033 ,0.0029 ,0.0026 ,0.0023 ,0.0021 ,0.0020 ,0.0019 ,0.0017 ,0.0016 ,0.0014 ,0.0012 ,0.0009 ,0.0006 ,0.0001 ,
                                                                0.0899,0.0816,0.0736,0.0660,0.0589,0.0524,0.0466,0.0415,0.0370,0.0331 ,0.0298 ,0.0269 ,0.0244 ,0.0222 ,0.0202 ,0.0184 ,0.0168 ,0.0154 ,0.0141 ,0.0129 ,0.0119 ,0.0110 ,0.0102 ,0.0096 ,0.0091 ,0.0086 ,0.0082 ,0.0078 ,0.0074 ,0.0070 ,0.0066 ,0.0061 ,0.0056 ,0.0050 ,0.0045 ,0.0039 ,0.0035 ,0.0031 ,0.0027 ,0.0025 ,0.0023 ,0.0021 ,0.0020 ,0.0018 ,0.0017 ,0.0015 ,0.0012 ,0.0009 ,0.0006 ,0.0001 ,
                                                                0.0916,0.0832,0.0752,0.0676,0.0604,0.0539,0.0480,0.0428,0.0383,0.0343 ,0.0309 ,0.0279 ,0.0253 ,0.0230 ,0.0210 ,0.0192 ,0.0175 ,0.0160 ,0.0147 ,0.0135 ,0.0124 ,0.0115 ,0.0107 ,0.0101 ,0.0095 ,0.0090 ,0.0086 ,0.0082 ,0.0078 ,0.0074 ,0.0069 ,0.0064 ,0.0058 ,0.0053 ,0.0047 ,0.0042 ,0.0037 ,0.0033 ,0.0029 ,0.0026 ,0.0024 ,0.0022 ,0.0021 ,0.0019 ,0.0017 ,0.0015 ,0.0013 ,0.0010 ,0.0006 ,0.0001 ,
                                                                0.0932,0.0849,0.0768,0.0692,0.0620,0.0554,0.0495,0.0442,0.0396,0.0355 ,0.0320 ,0.0289 ,0.0263 ,0.0239 ,0.0218 ,0.0199 ,0.0182 ,0.0167 ,0.0153 ,0.0141 ,0.0130 ,0.0120 ,0.0112 ,0.0106 ,0.0100 ,0.0095 ,0.0091 ,0.0086 ,0.0082 ,0.0078 ,0.0073 ,0.0067 ,0.0061 ,0.0055 ,0.0049 ,0.0044 ,0.0039 ,0.0034 ,0.0031 ,0.0028 ,0.0025 ,0.0023 ,0.0022 ,0.0020 ,0.0018 ,0.0016 ,0.0013 ,0.0010 ,0.0006 ,0.0001 ,
                                                                0.0948,0.0865,0.0785,0.0708,0.0636,0.0569,0.0509,0.0456,0.0408,0.0367 ,0.0331 ,0.0300 ,0.0272 ,0.0248 ,0.0227 ,0.0207 ,0.0189 ,0.0174 ,0.0159 ,0.0147 ,0.0136 ,0.0126 ,0.0118 ,0.0111 ,0.0105 ,0.0100 ,0.0095 ,0.0090 ,0.0086 ,0.0081 ,0.0076 ,0.0070 ,0.0064 ,0.0058 ,0.0051 ,0.0046 ,0.0041 ,0.0036 ,0.0032 ,0.0029 ,0.0027 ,0.0025 ,0.0023 ,0.0021 ,0.0019 ,0.0016 ,0.0014 ,0.0010 ,0.0006 ,0.0001 ,
                                                                0.0964,0.0882,0.0801,0.0724,0.0652,0.0585,0.0524,0.0470,0.0422,0.0379 ,0.0343 ,0.0311 ,0.0282 ,0.0257 ,0.0235 ,0.0215 ,0.0197 ,0.0181 ,0.0166 ,0.0153 ,0.0141 ,0.0132 ,0.0123 ,0.0116 ,0.0110 ,0.0104 ,0.0099 ,0.0095 ,0.0090 ,0.0085 ,0.0079 ,0.0073 ,0.0067 ,0.0060 ,0.0054 ,0.0048 ,0.0042 ,0.0038 ,0.0034 ,0.0031 ,0.0028 ,0.0026 ,0.0024 ,0.0022 ,0.0019 ,0.0017 ,0.0014 ,0.0010 ,0.0006 ,0.0000 ,
                                                                0.0980,0.0898,0.0817,0.0740,0.0667,0.0600,0.0539,0.0484,0.0435,0.0392 ,0.0354 ,0.0321 ,0.0292 ,0.0267 ,0.0244 ,0.0223 ,0.0204 ,0.0188 ,0.0173 ,0.0159 ,0.0147 ,0.0137 ,0.0128 ,0.0121 ,0.0115 ,0.0109 ,0.0104 ,0.0099 ,0.0094 ,0.0088 ,0.0083 ,0.0076 ,0.0069 ,0.0063 ,0.0056 ,0.0050 ,0.0044 ,0.0039 ,0.0035 ,0.0032 ,0.0029 ,0.0027 ,0.0025 ,0.0022 ,0.0020 ,0.0017 ,0.0014 ,0.0010 ,0.0006 ,0.0000 ,
                                                                0.0996,0.0914,0.0834,0.0756,0.0683,0.0616,0.0554,0.0498,0.0448,0.0404 ,0.0366 ,0.0332 ,0.0303 ,0.0276 ,0.0253 ,0.0231 ,0.0212 ,0.0195 ,0.0179 ,0.0166 ,0.0154 ,0.0143 ,0.0134 ,0.0126 ,0.0120 ,0.0114 ,0.0108 ,0.0103 ,0.0098 ,0.0092 ,0.0086 ,0.0079 ,0.0072 ,0.0065 ,0.0058 ,0.0052 ,0.0046 ,0.0041 ,0.0037 ,0.0033 ,0.0030 ,0.0028 ,0.0025 ,0.0023 ,0.0021 ,0.0018 ,0.0014 ,0.0011 ,0.0006 ,-0.0000,
                                                                0.1012,0.0930,0.0850,0.0773,0.0699,0.0631,0.0568,0.0512,0.0462,0.0417 ,0.0378 ,0.0343 ,0.0313 ,0.0286 ,0.0262 ,0.0240 ,0.0220 ,0.0202 ,0.0186 ,0.0172 ,0.0160 ,0.0149 ,0.0140 ,0.0132 ,0.0125 ,0.0119 ,0.0113 ,0.0108 ,0.0102 ,0.0096 ,0.0089 ,0.0082 ,0.0075 ,0.0067 ,0.0060 ,0.0054 ,0.0048 ,0.0043 ,0.0038 ,0.0035 ,0.0032 ,0.0029 ,0.0026 ,0.0024 ,0.0021 ,0.0018 ,0.0015 ,0.0011 ,0.0006 ,-0.0000,
                                                                0.1028,0.0946,0.0866,0.0789,0.0715,0.0647,0.0583,0.0526,0.0475,0.0430 ,0.0390 ,0.0355 ,0.0324 ,0.0296 ,0.0271 ,0.0248 ,0.0228 ,0.0210 ,0.0194 ,0.0179 ,0.0166 ,0.0155 ,0.0146 ,0.0137 ,0.0130 ,0.0124 ,0.0118 ,0.0112 ,0.0106 ,0.0099 ,0.0093 ,0.0085 ,0.0077 ,0.0070 ,0.0062 ,0.0056 ,0.0050 ,0.0044 ,0.0040 ,0.0036 ,0.0033 ,0.0030 ,0.0027 ,0.0025 ,0.0022 ,0.0019 ,0.0015 ,0.0011 ,0.0006 ,-0.0001,
                                                                0.1044,0.0963,0.0883,0.0805,0.0731,0.0662,0.0599,0.0541,0.0489,0.0443 ,0.0402 ,0.0366 ,0.0334 ,0.0306 ,0.0280 ,0.0257 ,0.0236 ,0.0218 ,0.0201 ,0.0186 ,0.0173 ,0.0162 ,0.0152 ,0.0143 ,0.0136 ,0.0129 ,0.0123 ,0.0116 ,0.0110 ,0.0103 ,0.0096 ,0.0088 ,0.0080 ,0.0072 ,0.0065 ,0.0058 ,0.0051 ,0.0046 ,0.0041 ,0.0037 ,0.0034 ,0.0031 ,0.0028 ,0.0025 ,0.0022 ,0.0019 ,0.0015 ,0.0011 ,0.0005 ,-0.0001,
                                                                0.1059,0.0979,0.0899,0.0821,0.0747,0.0678,0.0614,0.0555,0.0503,0.0456 ,0.0415 ,0.0378 ,0.0345 ,0.0316 ,0.0290 ,0.0266 ,0.0245 ,0.0226 ,0.0208 ,0.0193 ,0.0180 ,0.0168 ,0.0158 ,0.0149 ,0.0141 ,0.0134 ,0.0128 ,0.0121 ,0.0114 ,0.0107 ,0.0099 ,0.0091 ,0.0083 ,0.0075 ,0.0067 ,0.0059 ,0.0053 ,0.0047 ,0.0042 ,0.0038 ,0.0035 ,0.0032 ,0.0029 ,0.0026 ,0.0023 ,0.0019 ,0.0015 ,0.0011 ,0.0005 ,-0.0001};
        
        return arrayCurvatureMin[y * arrayWidth + x];
    }
    case 11:
    {
        static const std::vector<double> arrayCurvatureMin   = {0.0326,0.0244,0.0180,0.0131,0.0096,0.0070,0.0051,0.0038,0.0029,0.0022,0.0017,0.0013,0.0010,0.0008,0.0007,0.0005,0.0004,0.0004,0.0003,0.0003,0.0002,0.0002,0.0002,0.0001,0.0001,0.0001,0.0001,0.0000,0.0000,0.0000 ,-0.0000,-0.0000,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0008,-0.0008,-0.0009,
                                                                0.0348,0.0263,0.0195,0.0144,0.0106,0.0078,0.0058,0.0043,0.0032,0.0025,0.0019,0.0015,0.0012,0.0009,0.0007,0.0006,0.0005,0.0004,0.0003,0.0003,0.0002,0.0002,0.0002,0.0001,0.0001,0.0001,0.0001,0.0000,0.0000,0.0000 ,-0.0000,-0.0001,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0007,-0.0007,-0.0008,-0.0009,-0.0009,
                                                                0.0370,0.0282,0.0212,0.0158,0.0117,0.0087,0.0065,0.0049,0.0037,0.0028,0.0022,0.0017,0.0013,0.0010,0.0008,0.0007,0.0005,0.0004,0.0004,0.0003,0.0003,0.0002,0.0002,0.0002,0.0001,0.0001,0.0001,0.0000,0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0007,-0.0007,-0.0008,-0.0009,-0.0010,
                                                                0.0392,0.0302,0.0229,0.0173,0.0130,0.0097,0.0073,0.0056,0.0042,0.0033,0.0025,0.0019,0.0015,0.0012,0.0010,0.0008,0.0006,0.0005,0.0004,0.0003,0.0003,0.0002,0.0002,0.0002,0.0001,0.0001,0.0001,0.0000,0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0009,-0.0010,
                                                                0.0415,0.0323,0.0247,0.0188,0.0143,0.0109,0.0083,0.0064,0.0049,0.0038,0.0029,0.0023,0.0018,0.0014,0.0011,0.0009,0.0007,0.0006,0.0005,0.0004,0.0003,0.0003,0.0002,0.0002,0.0001,0.0001,0.0001,0.0000,0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0008,-0.0008,-0.0009,-0.0010,
                                                                0.0437,0.0343,0.0265,0.0204,0.0157,0.0121,0.0094,0.0073,0.0056,0.0044,0.0034,0.0027,0.0021,0.0017,0.0013,0.0010,0.0008,0.0007,0.0005,0.0004,0.0004,0.0003,0.0002,0.0002,0.0001,0.0001,0.0001,0.0000,0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0008,-0.0008,-0.0009,-0.0010,
                                                                0.0459,0.0363,0.0283,0.0220,0.0171,0.0134,0.0105,0.0083,0.0065,0.0051,0.0040,0.0032,0.0025,0.0020,0.0016,0.0012,0.0010,0.0008,0.0006,0.0005,0.0004,0.0003,0.0003,0.0002,0.0002,0.0001,0.0001,0.0000,0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0008,-0.0008,-0.0009,-0.0010,
                                                                0.0480,0.0383,0.0301,0.0236,0.0186,0.0147,0.0117,0.0093,0.0074,0.0059,0.0047,0.0038,0.0030,0.0024,0.0019,0.0015,0.0012,0.0009,0.0007,0.0006,0.0005,0.0004,0.0003,0.0002,0.0002,0.0001,0.0001,0.0001,0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0007,-0.0007,-0.0008,-0.0009,-0.0009,-0.0010,
                                                                0.0501,0.0402,0.0319,0.0252,0.0200,0.0159,0.0128,0.0104,0.0084,0.0068,0.0055,0.0045,0.0036,0.0029,0.0023,0.0018,0.0014,0.0011,0.0009,0.0007,0.0006,0.0004,0.0003,0.0003,0.0002,0.0002,0.0001,0.0001,0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0004,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0009,-0.0009,-0.0010,
                                                                0.0522,0.0422,0.0336,0.0267,0.0214,0.0172,0.0140,0.0114,0.0094,0.0078,0.0064,0.0052,0.0043,0.0035,0.0028,0.0022,0.0018,0.0014,0.0011,0.0009,0.0007,0.0005,0.0004,0.0003,0.0002,0.0002,0.0001,0.0001,0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0008,-0.0009,-0.0010,
                                                                0.0543,0.0441,0.0354,0.0283,0.0227,0.0184,0.0151,0.0125,0.0104,0.0087,0.0073,0.0061,0.0050,0.0041,0.0034,0.0027,0.0022,0.0017,0.0014,0.0011,0.0008,0.0007,0.0005,0.0004,0.0003,0.0002,0.0002,0.0001,0.0000,-0.0000,-0.0000,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0008,-0.0009,-0.0010,
                                                                0.0563,0.0460,0.0371,0.0298,0.0241,0.0196,0.0162,0.0135,0.0114,0.0096,0.0082,0.0069,0.0058,0.0049,0.0040,0.0033,0.0027,0.0021,0.0017,0.0013,0.0011,0.0008,0.0006,0.0005,0.0004,0.0003,0.0002,0.0001,0.0001,0.0000 ,-0.0000,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0008,-0.0009,-0.0010,
                                                                0.0584,0.0480,0.0389,0.0314,0.0255,0.0208,0.0173,0.0145,0.0123,0.0105,0.0090,0.0077,0.0066,0.0056,0.0047,0.0039,0.0032,0.0026,0.0021,0.0017,0.0013,0.0010,0.0008,0.0006,0.0005,0.0003,0.0002,0.0002,0.0001,0.0000 ,-0.0000,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0008,-0.0009,-0.0009,
                                                                0.0604,0.0499,0.0407,0.0330,0.0268,0.0220,0.0183,0.0155,0.0132,0.0114,0.0098,0.0085,0.0074,0.0064,0.0055,0.0047,0.0039,0.0032,0.0026,0.0021,0.0017,0.0013,0.0010,0.0008,0.0006,0.0005,0.0003,0.0002,0.0001,0.0001 ,0.0000 ,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0008,-0.0009,-0.0009,
                                                                0.0624,0.0518,0.0424,0.0345,0.0282,0.0233,0.0194,0.0164,0.0141,0.0122,0.0106,0.0093,0.0082,0.0072,0.0062,0.0054,0.0046,0.0039,0.0032,0.0026,0.0021,0.0017,0.0013,0.0010,0.0008,0.0006,0.0004,0.0003,0.0002,0.0001 ,0.0000 ,-0.0000,-0.0001,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0008,-0.0009,
                                                                0.0644,0.0538,0.0442,0.0361,0.0296,0.0245,0.0205,0.0174,0.0149,0.0130,0.0114,0.0100,0.0089,0.0079,0.0069,0.0061,0.0053,0.0046,0.0039,0.0032,0.0026,0.0021,0.0017,0.0013,0.0010,0.0008,0.0006,0.0004,0.0003,0.0002 ,0.0001 ,0.0000 ,-0.0001,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0008,-0.0009,
                                                                0.0663,0.0557,0.0460,0.0377,0.0311,0.0258,0.0216,0.0184,0.0158,0.0137,0.0121,0.0107,0.0095,0.0085,0.0076,0.0068,0.0060,0.0052,0.0045,0.0039,0.0032,0.0027,0.0022,0.0017,0.0013,0.0010,0.0008,0.0006,0.0004,0.0003 ,0.0001 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0009,
                                                                0.0682,0.0575,0.0478,0.0394,0.0325,0.0270,0.0227,0.0194,0.0167,0.0145,0.0128,0.0114,0.0102,0.0091,0.0082,0.0074,0.0066,0.0059,0.0052,0.0045,0.0039,0.0033,0.0027,0.0022,0.0017,0.0014,0.0010,0.0008,0.0006,0.0004 ,0.0002 ,0.0001 ,0.0000 ,-0.0001,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,-0.0008,
                                                                0.0701,0.0594,0.0495,0.0410,0.0340,0.0283,0.0239,0.0204,0.0176,0.0153,0.0135,0.0120,0.0108,0.0097,0.0088,0.0080,0.0072,0.0065,0.0058,0.0052,0.0045,0.0039,0.0033,0.0027,0.0022,0.0018,0.0014,0.0011,0.0008,0.0006 ,0.0004 ,0.0002 ,0.0001 ,0.0000 ,-0.0001,-0.0001,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0004,-0.0005,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,
                                                                0.0719,0.0613,0.0513,0.0426,0.0354,0.0297,0.0251,0.0214,0.0185,0.0162,0.0143,0.0127,0.0114,0.0102,0.0093,0.0085,0.0077,0.0071,0.0064,0.0058,0.0051,0.0045,0.0039,0.0033,0.0028,0.0023,0.0018,0.0014,0.0011,0.0008 ,0.0006 ,0.0004 ,0.0002 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0007,-0.0008,
                                                                0.0737,0.0631,0.0531,0.0442,0.0369,0.0310,0.0263,0.0225,0.0195,0.0170,0.0150,0.0134,0.0120,0.0108,0.0098,0.0090,0.0082,0.0075,0.0069,0.0063,0.0057,0.0051,0.0045,0.0039,0.0034,0.0028,0.0023,0.0018,0.0014,0.0011 ,0.0008 ,0.0006 ,0.0004 ,0.0002 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0008,
                                                                0.0754,0.0649,0.0548,0.0459,0.0384,0.0323,0.0275,0.0236,0.0205,0.0179,0.0158,0.0141,0.0126,0.0114,0.0103,0.0094,0.0087,0.0080,0.0074,0.0068,0.0063,0.0057,0.0051,0.0046,0.0040,0.0034,0.0029,0.0023,0.0019,0.0015 ,0.0011 ,0.0008 ,0.0006 ,0.0004 ,0.0002 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0007,-0.0008,
                                                                0.0771,0.0666,0.0565,0.0475,0.0399,0.0337,0.0287,0.0247,0.0215,0.0188,0.0166,0.0148,0.0132,0.0119,0.0108,0.0099,0.0091,0.0084,0.0078,0.0072,0.0067,0.0062,0.0057,0.0051,0.0046,0.0040,0.0035,0.0029,0.0024,0.0019 ,0.0015 ,0.0012 ,0.0009 ,0.0006 ,0.0004 ,0.0002 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0002,-0.0003,-0.0003,-0.0004,-0.0004,-0.0005,-0.0005,-0.0006,-0.0007,-0.0008,
                                                                0.0788,0.0684,0.0582,0.0491,0.0413,0.0350,0.0299,0.0259,0.0225,0.0198,0.0175,0.0156,0.0139,0.0125,0.0113,0.0104,0.0095,0.0088,0.0082,0.0076,0.0071,0.0066,0.0061,0.0056,0.0051,0.0046,0.0041,0.0035,0.0030,0.0025 ,0.0020 ,0.0016 ,0.0012 ,0.0009 ,0.0006 ,0.0004 ,0.0003 ,0.0001 ,0.0000 ,-0.0001,-0.0001,-0.0002,-0.0002,-0.0003,-0.0004,-0.0004,-0.0005,-0.0006,-0.0007,-0.0008,
                                                                0.0804,0.0700,0.0599,0.0506,0.0428,0.0363,0.0312,0.0270,0.0236,0.0208,0.0184,0.0163,0.0146,0.0131,0.0119,0.0108,0.0099,0.0092,0.0085,0.0080,0.0075,0.0070,0.0065,0.0061,0.0056,0.0052,0.0046,0.0041,0.0036,0.0031 ,0.0025 ,0.0021 ,0.0016 ,0.0013 ,0.0010 ,0.0007 ,0.0005 ,0.0003 ,0.0002 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0002,-0.0003,-0.0004,-0.0005,-0.0006,-0.0007,-0.0008,
                                                                0.0819,0.0717,0.0615,0.0522,0.0442,0.0377,0.0324,0.0281,0.0247,0.0218,0.0193,0.0172,0.0154,0.0138,0.0125,0.0113,0.0104,0.0096,0.0089,0.0083,0.0078,0.0073,0.0069,0.0065,0.0061,0.0056,0.0052,0.0047,0.0042,0.0037 ,0.0031 ,0.0026 ,0.0022 ,0.0017 ,0.0013 ,0.0010 ,0.0008 ,0.0005 ,0.0004 ,0.0002 ,0.0001 ,-0.0000,-0.0001,-0.0002,-0.0003,-0.0004,-0.0004,-0.0006,-0.0007,-0.0008,
                                                                0.0834,0.0733,0.0631,0.0537,0.0457,0.0390,0.0336,0.0293,0.0257,0.0228,0.0202,0.0180,0.0161,0.0145,0.0131,0.0119,0.0108,0.0100,0.0092,0.0086,0.0081,0.0076,0.0072,0.0068,0.0064,0.0060,0.0056,0.0052,0.0047,0.0043 ,0.0037 ,0.0032 ,0.0027 ,0.0023 ,0.0018 ,0.0014 ,0.0011 ,0.0008 ,0.0006 ,0.0004 ,0.0003 ,0.0001 ,0.0000 ,-0.0001,-0.0002,-0.0003,-0.0004,-0.0005,-0.0007,-0.0008,
                                                                0.0849,0.0749,0.0647,0.0553,0.0471,0.0403,0.0348,0.0304,0.0268,0.0238,0.0212,0.0189,0.0170,0.0152,0.0137,0.0124,0.0113,0.0104,0.0096,0.0089,0.0084,0.0079,0.0075,0.0071,0.0067,0.0064,0.0060,0.0056,0.0052,0.0048 ,0.0043 ,0.0038 ,0.0033 ,0.0028 ,0.0024 ,0.0019 ,0.0015 ,0.0012 ,0.0009 ,0.0007 ,0.0005 ,0.0003 ,0.0001 ,0.0000 ,-0.0001,-0.0003,-0.0004,-0.0005,-0.0007,-0.0008,
                                                                0.0863,0.0764,0.0663,0.0568,0.0485,0.0416,0.0360,0.0315,0.0279,0.0248,0.0222,0.0199,0.0178,0.0160,0.0144,0.0131,0.0119,0.0109,0.0100,0.0093,0.0087,0.0082,0.0077,0.0073,0.0070,0.0067,0.0063,0.0060,0.0057,0.0053 ,0.0049 ,0.0044 ,0.0040 ,0.0035 ,0.0030 ,0.0025 ,0.0021 ,0.0017 ,0.0013 ,0.0010 ,0.0007 ,0.0005 ,0.0003 ,0.0001 ,-0.0000,-0.0002,-0.0004,-0.0005,-0.0007,-0.0008,
                                                                0.0877,0.0779,0.0678,0.0583,0.0498,0.0428,0.0372,0.0326,0.0289,0.0258,0.0232,0.0208,0.0187,0.0168,0.0152,0.0137,0.0124,0.0113,0.0104,0.0096,0.0090,0.0084,0.0080,0.0076,0.0072,0.0069,0.0066,0.0063,0.0060,0.0057 ,0.0053 ,0.0050 ,0.0045 ,0.0041 ,0.0036 ,0.0031 ,0.0026 ,0.0022 ,0.0018 ,0.0014 ,0.0011 ,0.0008 ,0.0005 ,0.0003 ,0.0001 ,-0.0001,-0.0003,-0.0005,-0.0007,-0.0009,
                                                                0.0890,0.0794,0.0693,0.0597,0.0512,0.0441,0.0383,0.0337,0.0300,0.0269,0.0241,0.0218,0.0196,0.0177,0.0159,0.0144,0.0130,0.0119,0.0109,0.0100,0.0093,0.0087,0.0082,0.0078,0.0074,0.0071,0.0068,0.0066,0.0063,0.0060 ,0.0057 ,0.0054 ,0.0051 ,0.0046 ,0.0042 ,0.0037 ,0.0033 ,0.0028 ,0.0023 ,0.0019 ,0.0015 ,0.0011 ,0.0008 ,0.0005 ,0.0002 ,0.0000 ,-0.0002,-0.0004,-0.0006,-0.0009,
                                                                0.0903,0.0808,0.0708,0.0612,0.0525,0.0453,0.0395,0.0348,0.0310,0.0279,0.0251,0.0227,0.0205,0.0185,0.0167,0.0151,0.0137,0.0124,0.0113,0.0104,0.0097,0.0090,0.0085,0.0080,0.0077,0.0073,0.0071,0.0068,0.0066,0.0063 ,0.0061 ,0.0058 ,0.0055 ,0.0052 ,0.0048 ,0.0043 ,0.0039 ,0.0034 ,0.0029 ,0.0024 ,0.0019 ,0.0015 ,0.0011 ,0.0008 ,0.0004 ,0.0002 ,-0.0001,-0.0004,-0.0006,-0.0008,
                                                                0.0915,0.0822,0.0723,0.0626,0.0539,0.0465,0.0406,0.0359,0.0320,0.0289,0.0261,0.0237,0.0215,0.0194,0.0176,0.0159,0.0144,0.0130,0.0119,0.0109,0.0100,0.0093,0.0088,0.0083,0.0079,0.0075,0.0073,0.0070,0.0068,0.0066 ,0.0064 ,0.0061 ,0.0059 ,0.0056 ,0.0053 ,0.0049 ,0.0045 ,0.0040 ,0.0035 ,0.0030 ,0.0025 ,0.0020 ,0.0015 ,0.0011 ,0.0007 ,0.0004 ,0.0000 ,-0.0003,-0.0005,-0.0008,
                                                                0.0927,0.0836,0.0737,0.0640,0.0552,0.0477,0.0417,0.0369,0.0330,0.0298,0.0271,0.0246,0.0224,0.0203,0.0184,0.0167,0.0151,0.0137,0.0124,0.0114,0.0105,0.0097,0.0091,0.0085,0.0081,0.0077,0.0075,0.0072,0.0070,0.0068 ,0.0066 ,0.0064 ,0.0062 ,0.0060 ,0.0057 ,0.0054 ,0.0050 ,0.0046 ,0.0041 ,0.0036 ,0.0030 ,0.0025 ,0.0020 ,0.0015 ,0.0010 ,0.0006 ,0.0002 ,-0.0001,-0.0004,-0.0007,
                                                                0.0939,0.0849,0.0751,0.0654,0.0565,0.0489,0.0427,0.0379,0.0340,0.0308,0.0280,0.0256,0.0233,0.0213,0.0193,0.0175,0.0158,0.0143,0.0130,0.0119,0.0109,0.0101,0.0094,0.0088,0.0083,0.0080,0.0076,0.0074,0.0072,0.0070 ,0.0068 ,0.0066 ,0.0065 ,0.0063 ,0.0060 ,0.0058 ,0.0054 ,0.0051 ,0.0046 ,0.0042 ,0.0036 ,0.0031 ,0.0025 ,0.0020 ,0.0015 ,0.0010 ,0.0005 ,0.0001 ,-0.0003,-0.0006,
                                                                0.0950,0.0862,0.0765,0.0667,0.0577,0.0500,0.0438,0.0388,0.0349,0.0317,0.0290,0.0265,0.0243,0.0222,0.0202,0.0184,0.0166,0.0151,0.0137,0.0125,0.0114,0.0105,0.0097,0.0091,0.0086,0.0082,0.0079,0.0076,0.0073,0.0072 ,0.0070 ,0.0068 ,0.0067 ,0.0065 ,0.0063 ,0.0061 ,0.0058 ,0.0055 ,0.0051 ,0.0047 ,0.0042 ,0.0037 ,0.0031 ,0.0025 ,0.0019 ,0.0014 ,0.0009 ,0.0004 ,-0.0001,-0.0005,
                                                                0.0961,0.0875,0.0779,0.0681,0.0590,0.0512,0.0448,0.0398,0.0358,0.0326,0.0299,0.0274,0.0252,0.0231,0.0211,0.0192,0.0175,0.0158,0.0144,0.0131,0.0119,0.0110,0.0101,0.0095,0.0089,0.0084,0.0081,0.0078,0.0075,0.0073 ,0.0072 ,0.0070 ,0.0069 ,0.0067 ,0.0065 ,0.0063 ,0.0061 ,0.0058 ,0.0055 ,0.0051 ,0.0047 ,0.0042 ,0.0036 ,0.0031 ,0.0025 ,0.0019 ,0.0013 ,0.0007 ,0.0002 ,-0.0003,
                                                                0.0972,0.0887,0.0792,0.0694,0.0603,0.0523,0.0458,0.0407,0.0367,0.0334,0.0307,0.0283,0.0261,0.0241,0.0221,0.0201,0.0183,0.0166,0.0151,0.0137,0.0125,0.0115,0.0106,0.0098,0.0092,0.0087,0.0083,0.0080,0.0077,0.0075 ,0.0073 ,0.0072 ,0.0070 ,0.0069 ,0.0067 ,0.0065 ,0.0063 ,0.0061 ,0.0058 ,0.0055 ,0.0051 ,0.0047 ,0.0042 ,0.0036 ,0.0030 ,0.0024 ,0.0018 ,0.0012 ,0.0006 ,0.0000 ,
                                                                0.0982,0.0899,0.0805,0.0708,0.0615,0.0535,0.0468,0.0416,0.0375,0.0343,0.0316,0.0292,0.0270,0.0250,0.0230,0.0211,0.0192,0.0175,0.0159,0.0144,0.0131,0.0120,0.0110,0.0102,0.0096,0.0090,0.0086,0.0082,0.0079,0.0077 ,0.0075 ,0.0073 ,0.0072 ,0.0070 ,0.0069 ,0.0067 ,0.0065 ,0.0063 ,0.0061 ,0.0058 ,0.0055 ,0.0051 ,0.0046 ,0.0041 ,0.0035 ,0.0029 ,0.0023 ,0.0016 ,0.0010 ,0.0004 ,
                                                                0.0992,0.0911,0.0818,0.0721,0.0628,0.0546,0.0478,0.0425,0.0383,0.0351,0.0324,0.0300,0.0279,0.0259,0.0239,0.0220,0.0201,0.0183,0.0167,0.0152,0.0138,0.0126,0.0116,0.0107,0.0099,0.0093,0.0088,0.0084,0.0081,0.0079 ,0.0076 ,0.0074 ,0.0073 ,0.0071 ,0.0070 ,0.0068 ,0.0067 ,0.0065 ,0.0063 ,0.0060 ,0.0057 ,0.0054 ,0.0050 ,0.0045 ,0.0040 ,0.0034 ,0.0028 ,0.0022 ,0.0015 ,0.0008 ,
                                                                0.1001,0.0922,0.0831,0.0734,0.0640,0.0557,0.0488,0.0434,0.0391,0.0358,0.0331,0.0308,0.0287,0.0268,0.0248,0.0229,0.0210,0.0192,0.0175,0.0159,0.0145,0.0132,0.0121,0.0112,0.0104,0.0097,0.0091,0.0087,0.0083,0.0080 ,0.0078 ,0.0076 ,0.0074 ,0.0072 ,0.0071 ,0.0069 ,0.0068 ,0.0066 ,0.0064 ,0.0062 ,0.0059 ,0.0056 ,0.0053 ,0.0049 ,0.0044 ,0.0039 ,0.0033 ,0.0027 ,0.0020 ,0.0013 ,
                                                                0.1010,0.0933,0.0843,0.0747,0.0653,0.0568,0.0498,0.0442,0.0399,0.0366,0.0339,0.0316,0.0296,0.0276,0.0257,0.0238,0.0220,0.0202,0.0184,0.0168,0.0153,0.0139,0.0127,0.0117,0.0108,0.0101,0.0095,0.0090,0.0086,0.0082 ,0.0080 ,0.0077 ,0.0075 ,0.0073 ,0.0072 ,0.0070 ,0.0068 ,0.0067 ,0.0065 ,0.0063 ,0.0061 ,0.0058 ,0.0055 ,0.0052 ,0.0048 ,0.0043 ,0.0038 ,0.0032 ,0.0025 ,0.0018 ,
                                                                0.1019,0.0944,0.0856,0.0760,0.0665,0.0579,0.0508,0.0451,0.0407,0.0373,0.0346,0.0323,0.0303,0.0284,0.0266,0.0248,0.0229,0.0211,0.0193,0.0176,0.0161,0.0146,0.0134,0.0122,0.0113,0.0105,0.0098,0.0093,0.0088,0.0084 ,0.0081 ,0.0079 ,0.0076 ,0.0074 ,0.0073 ,0.0071 ,0.0069 ,0.0067 ,0.0066 ,0.0064 ,0.0062 ,0.0059 ,0.0057 ,0.0054 ,0.0050 ,0.0046 ,0.0041 ,0.0036 ,0.0030 ,0.0023 ,
                                                                0.1027,0.0955,0.0868,0.0773,0.0677,0.0591,0.0517,0.0459,0.0414,0.0380,0.0353,0.0330,0.0311,0.0292,0.0274,0.0256,0.0238,0.0220,0.0202,0.0185,0.0169,0.0154,0.0140,0.0129,0.0118,0.0110,0.0102,0.0096,0.0091,0.0087 ,0.0083 ,0.0080 ,0.0078 ,0.0075 ,0.0073 ,0.0072 ,0.0070 ,0.0068 ,0.0066 ,0.0064 ,0.0062 ,0.0060 ,0.0058 ,0.0055 ,0.0052 ,0.0048 ,0.0044 ,0.0039 ,0.0034 ,0.0028 ,
                                                                0.1036,0.0965,0.0880,0.0785,0.0690,0.0602,0.0527,0.0467,0.0421,0.0386,0.0359,0.0337,0.0318,0.0300,0.0283,0.0265,0.0248,0.0230,0.0212,0.0194,0.0178,0.0162,0.0148,0.0135,0.0124,0.0115,0.0106,0.0100,0.0094,0.0089 ,0.0085 ,0.0082 ,0.0079 ,0.0076 ,0.0074 ,0.0072 ,0.0070 ,0.0068 ,0.0067 ,0.0065 ,0.0063 ,0.0061 ,0.0058 ,0.0056 ,0.0053 ,0.0050 ,0.0046 ,0.0042 ,0.0037 ,0.0031 ,
                                                                0.1044,0.0975,0.0891,0.0798,0.0702,0.0613,0.0537,0.0476,0.0428,0.0393,0.0365,0.0343,0.0324,0.0307,0.0290,0.0274,0.0257,0.0239,0.0221,0.0204,0.0187,0.0170,0.0156,0.0142,0.0130,0.0120,0.0111,0.0103,0.0097,0.0092 ,0.0087 ,0.0083 ,0.0080 ,0.0077 ,0.0075 ,0.0073 ,0.0071 ,0.0069 ,0.0067 ,0.0065 ,0.0063 ,0.0061 ,0.0059 ,0.0056 ,0.0054 ,0.0051 ,0.0048 ,0.0044 ,0.0039 ,0.0034 ,
                                                                0.1051,0.0984,0.0903,0.0810,0.0714,0.0624,0.0547,0.0484,0.0436,0.0399,0.0371,0.0349,0.0331,0.0314,0.0298,0.0282,0.0265,0.0248,0.0231,0.0213,0.0196,0.0179,0.0164,0.0149,0.0137,0.0125,0.0116,0.0107,0.0100,0.0094 ,0.0089 ,0.0085 ,0.0082 ,0.0078 ,0.0076 ,0.0073 ,0.0071 ,0.0069 ,0.0067 ,0.0065 ,0.0063 ,0.0061 ,0.0059 ,0.0057 ,0.0054 ,0.0052 ,0.0048 ,0.0045 ,0.0041 ,0.0037 ,
                                                                0.1058,0.0994,0.0914,0.0822,0.0727,0.0636,0.0557,0.0492,0.0443,0.0405,0.0377,0.0355,0.0337,0.0321,0.0305,0.0290,0.0274,0.0257,0.0240,0.0223,0.0205,0.0188,0.0172,0.0157,0.0144,0.0131,0.0121,0.0112,0.0104,0.0097 ,0.0092 ,0.0087 ,0.0083 ,0.0080 ,0.0077 ,0.0074 ,0.0072 ,0.0069 ,0.0067 ,0.0065 ,0.0063 ,0.0061 ,0.0059 ,0.0057 ,0.0054 ,0.0052 ,0.0049 ,0.0046 ,0.0042 ,0.0038 ,
                                                                0.1065,0.1003,0.0925,0.0834,0.0739,0.0647,0.0567,0.0501,0.0450,0.0411,0.0383,0.0361,0.0343,0.0327,0.0312,0.0297,0.0282,0.0266,0.0249,0.0232,0.0215,0.0198,0.0181,0.0165,0.0151,0.0138,0.0126,0.0116,0.0108,0.0101 ,0.0094 ,0.0089 ,0.0085 ,0.0081 ,0.0078 ,0.0075 ,0.0072 ,0.0070 ,0.0068 ,0.0065 ,0.0063 ,0.0061 ,0.0059 ,0.0057 ,0.0054 ,0.0052 ,0.0049 ,0.0046 ,0.0043 ,0.0039 ,
                                                                0.1072,0.1012,0.0935,0.0846,0.0751,0.0659,0.0577,0.0509,0.0457,0.0417,0.0388,0.0366,0.0348,0.0333,0.0319,0.0304,0.0290,0.0275,0.0259,0.0242,0.0224,0.0207,0.0190,0.0174,0.0158,0.0145,0.0132,0.0122,0.0112,0.0104 ,0.0097 ,0.0091 ,0.0087 ,0.0082 ,0.0079 ,0.0076 ,0.0073 ,0.0070 ,0.0068 ,0.0066 ,0.0063 ,0.0061 ,0.0059 ,0.0057 ,0.0054 ,0.0052 ,0.0049 ,0.0046 ,0.0043 ,0.0039};
        
        return arrayCurvatureMin[y * arrayWidth + x];
    }
    }
    
    return 0;
}
