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

void drawAOI(cv::Mat& I, const AOIProperties &mAOI, const cv::Vec3b &col)
{
    int imgWdth = I.cols;
    int imgHght = I.rows;

    int AOIXPos = mAOI.xPos;
    int AOIYPos = mAOI.yPos;
    int AOIWdth = mAOI.wdth;
    int AOIHght = mAOI.hght;

    if (AOIXPos < 0) { AOIXPos = 0; }
    if (AOIYPos < 0) { AOIYPos = 0; }
    if (AOIXPos + AOIWdth >= imgWdth) { AOIWdth = imgWdth - AOIXPos - 1; }
    if (AOIYPos + AOIHght >= imgHght) { AOIHght = imgHght - AOIYPos - 1; }

    for (int x = AOIXPos; x < AOIXPos + AOIWdth; x++) // left to right
    {
        I.at<cv::Vec3b>(AOIYPos,           x) = col; // top
        I.at<cv::Vec3b>(AOIYPos + AOIHght, x) = col; // bottom
    }

    for (int y = AOIYPos; y < AOIYPos + AOIHght; y++) // top to bottom
    {
        I.at<cv::Vec3b>(y, AOIXPos)           = col; // left
        I.at<cv::Vec3b>(y, AOIXPos + AOIWdth) = col; // right
    }
}

void drawEdges(cv::Mat& I, const AOIProperties& mAOI, const cv::Vec3b& col, const std::vector<int>& edgeIndices)
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

void drawOutline(cv::Mat& I, const AOIProperties& mAOI, const cv::Vec3b& colour_1, const cv::Vec3b& colour_2, const cv::Vec3b&  colour_3, const std::vector<edgeProperties>& vEdgePropertiesAll)
{
    int numEdges = vEdgePropertiesAll.size();

    for (int iEdge = 0; iEdge < numEdges; iEdge++)
    {
        int edgeSize = vEdgePropertiesAll[iEdge].pointIndices.size();

        cv::Vec3b colour;

        if      (vEdgePropertiesAll[iEdge].tag == 2) { colour = colour_1; }
        else if (vEdgePropertiesAll[iEdge].tag == 1) { colour = colour_2; }
        else                                         { colour = colour_3; }

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

void drawEllipse(cv::Mat& I, const AOIProperties& mAOI, const cv::Vec3b& col, const std::vector<double>& c)
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

void drawCross(cv::Mat& I, double crossX, double crossY, int ellipseDrawCrossSize, const cv::Vec3b& col)
{
    int wdth = I.cols;
    int hght = I.rows;

    if (crossX >= 0 && crossX < wdth && crossY >= 0 && crossY < hght)
    {
        for (int x = crossX - ellipseDrawCrossSize; x <= crossX + ellipseDrawCrossSize; x++)
        {
            int X = x;
            int Y = crossY;

            if (X >= 0 && X < wdth)
            {
                I.at<cv::Vec3b>(Y, X) = col;
            }
        }

        for (int y = crossY - ellipseDrawCrossSize; y <= crossY + ellipseDrawCrossSize; y++)
        {
            int X = crossX;
            int Y = y;

            if (Y >= 0 && Y < hght)
            {
                I.at<cv::Vec3b>(Y, X) = col;
            }
        }
    }
}

void drawAll(cv::Mat &I, const drawVariables &mDrawVariables)
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
            drawAOI(I, mDrawVariables.outerAOI, blue);
            drawAOI(I, mDrawVariables.innerAOI, blue);
            drawAOI(I, mDrawVariables.glintAOI, blue);
        }

        if (Parameters::drawFlags.edge)
        {
            drawEdges  (I, mDrawVariables.outerAOI, red, mDrawVariables.cannyEdgeIndices);
            drawOutline(I, mDrawVariables.outerAOI, green, yellow, orange, mDrawVariables.edgePropertiesAll);
        }

        drawCross(I, mDrawVariables.predictedXPos, mDrawVariables.predictedYPos, Parameters::ellipseDrawCrossSize, cyan);

        if (Parameters::drawFlags.elps)
        {
            if (mDrawVariables.DETECTED)
            {
                drawEllipse(I, mDrawVariables.outerAOI, white, mDrawVariables.ellipseCoefficients);
                drawCross(I, mDrawVariables.exactXPos, mDrawVariables.exactYPos, Parameters::ellipseDrawCrossSize, white);
            }
        }

    }
}
