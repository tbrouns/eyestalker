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

// draw functions

void drawHaarDetector(cv::Mat& I, int x0, int y0, int haarWidth, cv::Vec3b col)
{
    unsigned short int wdth = I.rows;
    unsigned short int hgth = I.cols;

    for (unsigned short int y = 0; y < wdth; y++)
    {
        for (unsigned short int x = 0; x < hgth; x++)
        {
            if (x == x0 || x == x0 + haarWidth)
            {
                if (y > y0 && y < y0 + haarWidth)
                {
                    I.at<cv::Vec3b>(y, x) = col;
                }
            }

            if (y == y0 || y == y0 + haarWidth)
            {
                if (x > x0 && x < x0 + haarWidth)
                {
                    I.at<cv::Vec3b>(y, x) = col;
                }
            }
        }
    }
}

void drawEdges(cv::Mat& I, const std::vector<char>& p, int x0, int y0, int haarWidth, const cv::Vec3b& col)
{
    int p_size = p.size();

    if (p_size > 0)
    {
        for (int y = 0; y < haarWidth; y++)
        {
            for (int x = 0; x < haarWidth; x++)
            {
                int i = y * haarWidth + x;

                int X = x + x0;
                int Y = y + y0;

                if (p[i] == 1)
                {
                    I.at<cv::Vec3b>(Y, X) = col;
                }
            }
        }
    }
}

void drawOutline(cv::Mat& I, const std::vector<edgeProperties>& vEdgePropertiesAll, int x0, int y0, int haarWidth, const cv::Vec3b& primaryColour, const cv::Vec3b& secondaryColour)
{
    int numEdges = vEdgePropertiesAll.size();

    for (int iEdge = 0; iEdge < numEdges; iEdge++)
    {
        int edgeSize = vEdgePropertiesAll[iEdge].size;

        cv::Vec3b colour;

        if (vEdgePropertiesAll[iEdge].flag) { colour = primaryColour; }
        else                                { colour = secondaryColour; }

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

void drawEllipse(cv::Mat& I, const std::vector<double>& c, int x0, int y0, int haarWidth, const cv::Vec3b& col)
{
    int c_size = c.size();

    double lineWidth = round(Parameters::ellipseDrawOutlineWidth * (I.cols + I.rows));

    if (c_size == 6)
    {
        for (int x = 0; x < haarWidth; x++)
        {
            for (int y = 0; y < haarWidth; y++)
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

void drawAll(cv::Mat &I, eyeProperties mEyeProperties)
{
    cv::Vec3b blue(255, 0, 0);
    cv::Vec3b green(0, 255, 0);
    cv::Vec3b red(0, 0, 255);

    cv::Vec3b cyan(255, 255, 0);
    cv::Vec3b yellow(0, 255, 255);

    cv::Vec3b white(255, 255, 255);


    if (Parameters::drawFlags.haar)
    {
        drawHaarDetector(I, mEyeProperties.m.offsetPupilHaarXPos, mEyeProperties.m.offsetPupilHaarYPos, mEyeProperties.m.offsetPupilHaarWdth, blue);
        drawHaarDetector(I, mEyeProperties.m.pupilHaarXPos, mEyeProperties.m.pupilHaarYPos, mEyeProperties.m.pupilHaarWdth, blue);
        drawHaarDetector(I, mEyeProperties.m.glintHaarXPos, mEyeProperties.m.glintHaarYPos, mEyeProperties.m.glintHaarWdth, blue);
    }

    if (Parameters::drawFlags.edge)
    {
        drawEdges(I, mEyeProperties.m.cannyEdges, mEyeProperties.m.offsetPupilHaarXPos, mEyeProperties.m.offsetPupilHaarYPos, mEyeProperties.m.offsetPupilHaarWdth, red);
        drawOutline(I, mEyeProperties.m.edgePropertiesAll, mEyeProperties.m.offsetPupilHaarXPos, mEyeProperties.m.offsetPupilHaarYPos, mEyeProperties.m.offsetPupilHaarWdth, green, yellow);
    }

    if (Parameters::drawFlags.elps)
    {
        if (mEyeProperties.v.pupilDetected)
        {
            drawEllipse(I, mEyeProperties.m.ellipseCoefficients, mEyeProperties.m.offsetPupilHaarXPos, mEyeProperties.m.offsetPupilHaarYPos, mEyeProperties.m.offsetPupilHaarWdth, white);
            drawEllipseCross(I, mEyeProperties.v.xPosExact, mEyeProperties.v.yPosExact, Parameters::ellipseDrawCrossSize, white);
        }
    }

    drawEllipseCross(I, mEyeProperties.v.xPosPredicted, mEyeProperties.v.yPosPredicted, Parameters::ellipseDrawCrossSize, cyan);
}

void drawROI(cv::Mat& I, int x0, int y0, int W, int H, cv::Vec3b col)
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
