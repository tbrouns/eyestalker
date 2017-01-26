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

void drawHaarDetector(cv::Mat& I, int x0, int y0, int haarWdth, int haarHght, cv::Vec3b col)
{
    int imgWdth = I.cols;
    int imgHght = I.rows;

    if (x0 < 0) { x0 = 0; }
    if (y0 < 0) { y0 = 0; }
    if (x0 + haarWdth >= imgWdth) { haarWdth = imgWdth - x0 - 1; }
    if (y0 + haarHght >= imgHght) { haarHght = imgHght - y0 - 1; }

    for (int x = x0; x < x0 + haarWdth; x++) // left to right
    {
        I.at<cv::Vec3b>(y0,            x) = col; // top
        I.at<cv::Vec3b>(y0 + haarHght, x) = col; // bottom
    }

    for (int y = y0; y < y0 + haarHght; y++) // top to bottom
    {
        I.at<cv::Vec3b>(y, x0)             = col; // left
        I.at<cv::Vec3b>(y, x0 + haarWdth) = col; // right
    }
}

void drawEdges(cv::Mat& I, const std::vector<int>& edgeIndices, int x0, int y0, int haarWidth, const cv::Vec3b& col)
{
    int numEdgePoints = edgeIndices.size();

    for (int iEdgePoint = 0; iEdgePoint < numEdgePoints; iEdgePoint++)
    {
        int edgePointIndex = edgeIndices[iEdgePoint];
        int x =  edgePointIndex % haarWidth;
        int y = (edgePointIndex - x) / haarWidth;

        int X = x + x0;
        int Y = y + y0;

        I.at<cv::Vec3b>(Y, X) = col;
    }
}

void drawOutline(cv::Mat& I, const std::vector<edgeProperties>& vEdgePropertiesAll, int x0, int y0, int haarWidth, const cv::Vec3b& primaryColour, const cv::Vec3b& secondaryColour, const cv::Vec3b&  tertiaryColour)
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

            int x = edgePointIndex % haarWidth;
            int y = (edgePointIndex - x) / haarWidth;

            int X = x + x0;
            int Y = y + y0;

            I.at<cv::Vec3b>(Y, X) = colour;
        }
    }
}

void drawEllipse(cv::Mat& I, const std::vector<double>& c, int x0, int y0, int haarWidth, int haarHeight, const cv::Vec3b& col)
{
    int c_size = c.size();

    double lineWidth = round(Parameters::ellipseDrawOutlineWidth * (I.cols + I.rows));

    if (c_size == 6)
    {
        for (int x = 0; x < haarWidth; x++)
        {
            for (int y = 0; y < haarHeight; y++)
            {
                double b = c[0] * x * x + c[1] * x * y + c[2] * y * y + c[3] * x + c[4] * y + c[5];

                if (b < lineWidth && b > -lineWidth)
                {
                    int X = x + x0;
                    int Y = y + y0;

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
    cv::Vec3b yellow(  0, 255, 255);
    cv::Vec3b white (255, 255, 255);

    if (I.cols > 0 && I.rows > 0)
    {
        if (Parameters::drawFlags.haar)
        {
            drawHaarDetector(I, mDetectionProperties.m.offsetPupilHaarXPos, mDetectionProperties.m.offsetPupilHaarYPos, mDetectionProperties.m.offsetPupilHaarWdth, mDetectionProperties.m.offsetPupilHaarHght, blue);
            drawHaarDetector(I, mDetectionProperties.m.pupilHaarXPos, mDetectionProperties.m.pupilHaarYPos, mDetectionProperties.m.pupilHaarWdth, mDetectionProperties.m.pupilHaarHght, blue);
            drawHaarDetector(I, mDetectionProperties.m.glintXPos, mDetectionProperties.m.glintYPos, mDetectionProperties.m.glintSize, mDetectionProperties.m.glintSize, blue);
        }

        if (Parameters::drawFlags.edge)
        {
            drawEdges(I, mDetectionProperties.m.cannyEdgeIndices, mDetectionProperties.m.offsetPupilHaarXPos, mDetectionProperties.m.offsetPupilHaarYPos, mDetectionProperties.m.offsetPupilHaarWdth, red);
            drawOutline(I, mDetectionProperties.m.edgePropertiesAll, mDetectionProperties.m.offsetPupilHaarXPos, mDetectionProperties.m.offsetPupilHaarYPos, mDetectionProperties.m.offsetPupilHaarWdth, cyan, green, yellow);
        }

        if (Parameters::drawFlags.elps)
        {
            if (mDetectionProperties.v.pupilDetected)
            {
                drawEllipse(I, mDetectionProperties.m.ellipseCoefficients, mDetectionProperties.m.offsetPupilHaarXPos, mDetectionProperties.m.offsetPupilHaarYPos, mDetectionProperties.m.offsetPupilHaarWdth, mDetectionProperties.m.offsetPupilHaarHght, white);
                drawEllipseCross(I, mDetectionProperties.v.xPosExact, mDetectionProperties.v.yPosExact, Parameters::ellipseDrawCrossSize, white);
            }
        }

        drawEllipseCross(I, mDetectionProperties.v.xPosPredicted, mDetectionProperties.v.yPosPredicted, Parameters::ellipseDrawCrossSize, cyan);
    }
}

void drawAOI(cv::Mat& I, int x0, int y0, int W, int H, cv::Vec3b col)
{
    int wdth = I.rows;
    int hgth = I.cols;

    for (int y = 0; y < wdth; y++)
    {
        for (int x = 0; x < hgth; x++)
        {
            if (x == x0 || x == x0 + W)
            {
                if (y > y0 && y < y0 + H)
                {
                    I.at<cv::Vec3b>(y, x) = col;
                }
            }

            if (y == y0 || y == y0 + H)
            {
                if (x > x0 && x < x0 + W)
                {
                    I.at<cv::Vec3b>(y, x) = col;
                }
            }
        }
    }
}
