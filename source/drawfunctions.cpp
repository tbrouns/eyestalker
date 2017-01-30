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

#include "headers/drawfunctions.h"

void drawAOI(cv::Mat& I, AOIProperties mAOI, cv::Vec3b col)
{
    int imgWdth = I.cols;
    int imgHght = I.rows;

    if (mAOI.xPos < 0) { mAOI.xPos = 0; }
    if (mAOI.yPos < 0) { mAOI.yPos = 0; }
    if (mAOI.xPos + mAOI.wdth >= imgWdth) { mAOI.wdth = imgWdth - mAOI.xPos - 1; }
    if (mAOI.yPos + mAOI.hght >= imgHght) { mAOI.hght = imgHght - mAOI.yPos - 1; }

    for (int x = mAOI.xPos; x < mAOI.xPos + mAOI.wdth; x++) // left to right
    {
        I.at<cv::Vec3b>(mAOI.yPos,             x) = col; // top
        I.at<cv::Vec3b>(mAOI.yPos + mAOI.hght, x) = col; // bottom
    }

    for (int y = mAOI.yPos; y < mAOI.yPos + mAOI.hght; y++) // top to bottom
    {
        I.at<cv::Vec3b>(y, mAOI.xPos)             = col; // left
        I.at<cv::Vec3b>(y, mAOI.xPos + mAOI.wdth) = col; // right
    }
}

void drawEdges(cv::Mat& I, const std::vector<int>& edgeIndices, AOIProperties mAOI, const cv::Vec3b& col)
{
    int numEdgePoints = edgeIndices.size();

    for (int iEdgePoint = 0; iEdgePoint < numEdgePoints; iEdgePoint++)
    {
        int edgePointIndex = edgeIndices[iEdgePoint];
        int x =  edgePointIndex % mAOI.wdth;
        int y = (edgePointIndex - x) / mAOI.wdth;

        int X = x + mAOI.xPos;
        int Y = y + mAOI.yPos;

        I.at<cv::Vec3b>(Y, X) = col;
    }
}

void drawOutline(cv::Mat& I, const std::vector<edgeProperties>& vEdgePropertiesAll, AOIProperties mAOI, const cv::Vec3b& primaryColour, const cv::Vec3b& secondaryColour, const cv::Vec3b&  tertiaryColour)
{
    int numEdges = vEdgePropertiesAll.size();

    for (int iEdge = 0; iEdge < numEdges; iEdge++)
    {
        int edgeSize = vEdgePropertiesAll[iEdge].size;

        cv::Vec3b colour;

        if      (vEdgePropertiesAll[iEdge].flag == 2) { colour = primaryColour;   }
        else if (vEdgePropertiesAll[iEdge].flag == 1) { colour = secondaryColour; }
        else                                          { colour = tertiaryColour;  }

        for (int iEdgePoint = 0; iEdgePoint < edgeSize; iEdgePoint++)
        {
            int edgePointIndex = vEdgePropertiesAll[iEdge].pointIndices[iEdgePoint];

            int x = edgePointIndex % mAOI.wdth;
            int y = (edgePointIndex - x) / mAOI.wdth;

            int X = x + mAOI.xPos;
            int Y = y + mAOI.yPos;

            I.at<cv::Vec3b>(Y, X) = colour;
        }
    }
}

void drawEllipse(cv::Mat& I, const std::vector<double>& c, AOIProperties mAOI, const cv::Vec3b& col)
{
    int c_size = c.size();

    double lineWidth = round(Parameters::ellipseDrawOutlineWidth * (I.cols + I.rows));

    if (c_size == 6)
    {
        for (int x = 0; x < mAOI.wdth; x++)
        {
            for (int y = 0; y < mAOI.hght; y++)
            {
                double b = c[0] * x * x + c[1] * x * y + c[2] * y * y + c[3] * x + c[4] * y + c[5];

                if (b < lineWidth && b > -lineWidth)
                {
                    int X = x + mAOI.xPos;
                    int Y = y + mAOI.yPos;

                    I.at<cv::Vec3b>(Y, X) = col;
                }
            }
        }
    }
}

void drawEllipseCross(cv::Mat& I, double cx, double cy, int ellipseDrawCrossSize, const cv::Vec3b& col)
{
    int wdth = I.cols;
    int hght = I.rows;

    for (int x = round(cx) - ellipseDrawCrossSize; x <= round(cx) + ellipseDrawCrossSize; x++)
    {
        int X = x;
        int Y = round(cy);

        if (X >= 0 && X < wdth)
        {
            I.at<cv::Vec3b>(Y, X) = col;
        }
    }

    for (int y = round(cy) - ellipseDrawCrossSize; y <= round(cy) + ellipseDrawCrossSize; y++)
    {
        int X = round(cx);
        int Y = y;

        if (Y >= 0 && Y < hght)
        {
            I.at<cv::Vec3b>(Y, X) = col;
        }
    }
}

void drawAll(cv::Mat &I, detectionProperties mDetectionProperties)
{
    cv::Vec3b blue  (255,   0,   0);
    cv::Vec3b green (  0, 255,   0);
    cv::Vec3b red   (  0,   0, 255);
    cv::Vec3b cyan  (255, 255,   0);
    cv::Vec3b orange(  0, 165, 255);
    cv::Vec3b yellow(  0, 255, 255);
    cv::Vec3b white (255, 255, 255);

    if (I.cols > 0 && I.rows > 0)
    {
        if (Parameters::drawFlags.haar)
        {
            drawAOI(I, mDetectionProperties.m.outerAOI, blue);
            drawAOI(I, mDetectionProperties.m.innerAOI, blue);
            drawAOI(I, mDetectionProperties.m.glintAOI, blue);
        }

        if (Parameters::drawFlags.edge)
        {
            drawEdges  (I, mDetectionProperties.m.cannyEdgeIndices,  mDetectionProperties.m.outerAOI, red);
            drawOutline(I, mDetectionProperties.m.edgePropertiesAll, mDetectionProperties.m.outerAOI, green, yellow, orange);
        }

        drawEllipseCross(I, mDetectionProperties.v.xPosPrediction, mDetectionProperties.v.yPosPrediction, Parameters::ellipseDrawCrossSize, cyan);

        if (Parameters::drawFlags.elps)
        {
            if (mDetectionProperties.v.pupilDetected)
            {
                drawEllipse(I, mDetectionProperties.m.ellipseCoefficients, mDetectionProperties.m.outerAOI, white);
                drawEllipseCross(I, mDetectionProperties.v.xPosExact, mDetectionProperties.v.yPosExact, Parameters::ellipseDrawCrossSize, white);
            }
        }

    }
}
