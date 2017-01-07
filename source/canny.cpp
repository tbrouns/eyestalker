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

#include "headers/canny.h"

void radialGradient(cv::Mat src, cv::Mat dstX, cv::Mat dstY, int aperture_size)
{
    int kernelPerimeter = 8; // perimeter length of kernel
    double a = 12.0;
    double b =  3.0;

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

    std::vector<double> wX(kernelPerimeter);
    wX[0] = 1.0;
    wX[1] = 0.5;
    wX[2] = 0.0;
    wX[3] = 0.5;
    wX[4] = 1.0;
    wX[5] = 0.5;
    wX[6] = 0.0;
    wX[7] = 0.5;

    std::vector<double> wY(kernelPerimeter);
    wY[0] = 0.0;
    wY[1] = 0.5;
    wY[2] = 1.0;
    wY[3] = 0.5;
    wY[4] = 0.0;
    wY[5] = 0.5;
    wY[6] = 1.0;
    wY[7] = 0.5;

    uchar *ptr = src.data;
    int width  = src.cols;
    int height = src.rows;

    int borderSize = (aperture_size - 1) / 2;

    int xCentre = 0.5 * width;
    int yCentre = 0.5 * height;

    for (int y = borderSize; y < height - borderSize; y++)
    {
        for (int x = borderSize; x < width - borderSize; x++)
        {
            double dx = x - xCentre;
            double dy = yCentre - y;

            double norm = sqrt(pow(dx,2) + pow(dy,2));
            dx = dx / norm;
            dy = dy / norm;

            double theta;
            if (dx != 0 && dy != 0)
            {
                theta = std::abs(atan(dy/dx));
                if      (dx > 0 && dy > 0) { theta = theta; }
                else if (dx > 0 && dy < 0) { theta = 2 * M_PI - theta; }
                else if (dx < 0 && dy > 0) { theta = M_PI - theta; }
                else if (dx < 0 && dy < 0) { theta = M_PI + theta; }
            }
            else if (dx == 0 && dy != 0)
            {
                if (dy > 0) { theta = 0.5 * M_PI; }
                else        { theta = 1.5 * M_PI; }
            }
            else if (dx != 0 && dy == 0)
            {
                if (dx > 0) { theta = 0;    }
                else        { theta = M_PI; }
            }
            else { theta = 0.5 * M_PI; }   // arbitrary choice

            double alpha = theta * kernelPerimeter / (2 * M_PI);
            double valX = 0;
            double valY = 0;

            for (int i = 0; i < kernelPerimeter; i++)
            {
                double dpos = std::abs(i - alpha);
                if (dpos >  kernelPerimeter * 0.5) { dpos = kernelPerimeter - dpos; }
                double dneg = 0.5 * kernelPerimeter - dpos;

                double kernelVal = a * (exp(- pow(dpos, 2) / b) - exp(- pow(dneg, 2) / b));

                // Do convolution
                int iPixel = width * (y + dY[i] * borderSize) + (x + dX[i] * borderSize);
                double  pixelIntensity = ptr[iPixel];
                valX += pixelIntensity * wX[i] * kernelVal;
                valY += pixelIntensity * wY[i] * kernelVal;
            }

            dstX.at<short>(y, x) = round(valX);
            dstY.at<short>(y, x) = round(valY);
        }
    }
}


/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                        Intel License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000, Intel Corporation, all rights reserved.
// Copyright (C) 2014, Itseez Inc., all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of Intel Corporation may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

cv::Mat cannyEdgeDetection(cv::Mat src, double low_thresh, double high_thresh, int aperture_size, bool L2gradient)
{
    const int type = src.type(), depth = CV_MAT_DEPTH(type), cn = CV_MAT_CN(type);

    if ((aperture_size & 1) == 0 || (aperture_size != -1 && (aperture_size < 3 || aperture_size > 7)))
        CV_Error(CV_StsBadFlag, "Aperture size should be odd");

    if (low_thresh > high_thresh)
        std::swap(low_thresh, high_thresh);

    int borderSize = (aperture_size - 1) / 2;

    cv::Mat dx(src.rows, src.cols, CV_16SC(cn));
    cv::Mat dy(src.rows, src.cols, CV_16SC(cn));
    radialGradient(src, dx, dy, aperture_size);

//    cv::Mat dx(src.rows, src.cols, CV_16SC(cn));
//    cv::Mat dy(src.rows, src.cols, CV_16SC(cn));
//    cv::Sobel(src, dx, CV_16S, 1, 0, aperture_size, 1, 0, cv::BORDER_REPLICATE);
//    cv::Sobel(src, dy, CV_16S, 0, 1, aperture_size, 1, 0, cv::BORDER_REPLICATE);

    cv::Rect  cropAOI(borderSize, borderSize, src.rows - 2 * borderSize, src.cols - 2 * borderSize);
    dx  =  dx(cropAOI);
    dy  =  dy(cropAOI);
    src = src(cropAOI);

    const cv::Size size = src.size();

    cv::Mat dst;
    CV_Assert( depth == CV_8U );
    dst.create(size, CV_8U);

    cv::imwrite( "G:/Media/Pictures/eyes/dx.jpg",  dx);
    cv::imwrite( "G:/Media/Pictures/eyes/dy.jpg",  dy);
//    cv::imwrite( "G:/Media/Pictures/eyes/dx2.jpg", dx2);
//    cv::imwrite( "G:/Media/Pictures/eyes/dy2.jpg", dy2);
    cv::imwrite( "G:/Media/Pictures/eyes/src.jpg", src);

    if (L2gradient)
    {
        low_thresh  = std::min(32767.0,  low_thresh);
        high_thresh = std::min(32767.0, high_thresh);

        if ( low_thresh > 0)  low_thresh *=  low_thresh;
        if (high_thresh > 0) high_thresh *= high_thresh;
    }
    int low  = cvFloor( low_thresh);
    int high = cvFloor(high_thresh);

    ptrdiff_t mapstep = src.cols + 2;
    cv::AutoBuffer<uchar> buffer((src.cols+2)*(src.rows+2) + cn * mapstep * 3 * sizeof(int));

    int* mag_buf[3];
    mag_buf[0] = (int*)(uchar*)buffer;
    mag_buf[1] = mag_buf[0] + mapstep*cn;
    mag_buf[2] = mag_buf[1] + mapstep*cn;
    memset(mag_buf[0], 0, /* cn* */mapstep*sizeof(int));

    uchar* map = (uchar*)(mag_buf[2] + mapstep*cn);
    memset(map, 1, mapstep);
    memset(map + mapstep*(src.rows + 1), 1, mapstep);

    int maxsize = std::max(1 << 10, src.cols * src.rows / 10);
    std::vector<uchar*> stack(maxsize);
    uchar **stack_top = &stack[0];
    uchar **stack_bottom = &stack[0];

    /* sector numbers
       (Top-Left Origin)

        1   2   3
         *  *  *
          * * *
        0*******0
          * * *
         *  *  *
        3   2   1
    */

#define CANNY_PUSH(d)    *(d) = uchar(2), *stack_top++ = (d)
#define CANNY_POP(d)      (d) = *--stack_top

    // calculate magnitude and angle of gradient, perform non-maxima suppression.
    // fill the map with one of the following values:
    //   0 - the pixel might belong to an edge
    //   1 - the pixel can not belong to an edge
    //   2 - the pixel does belong to an edge
    for (int i = 0; i <= src.rows; i++)
    {
        int* _norm = mag_buf[(i > 0) + 1] + 1;
        if (i < src.rows)
        {
            short* _dx = dx.ptr<short>(i);
            short* _dy = dy.ptr<short>(i);

            if (!L2gradient)
            {
                int j = 0, width = src.cols * cn;
                for ( ; j < width; ++j)
                    _norm[j] = std::abs(int(_dx[j])) + std::abs(int(_dy[j]));
            }
            else
            {
                int j = 0, width = src.cols * cn;
                for ( ; j < width; ++j)
                    _norm[j] = int(_dx[j])*_dx[j] + int(_dy[j])*_dy[j];
            }

            if (cn > 1)
            {
                for(int j = 0, jn = 0; j < src.cols; ++j, jn += cn)
                {
                    int maxIdx = jn;
                    for(int k = 1; k < cn; ++k)
                        if(_norm[jn + k] > _norm[maxIdx]) maxIdx = jn + k;
                    _norm[j] = _norm[maxIdx];
                    _dx[j] = _dx[maxIdx];
                    _dy[j] = _dy[maxIdx];
                }
            }
            _norm[-1] = _norm[src.cols] = 0;
        }
        else
            memset(_norm-1, 0, /* cn* */mapstep*sizeof(int));

        // at the very beginning we do not have a complete ring
        // buffer of 3 magnitude rows for non-maxima suppression
        if (i == 0)
            continue;

        uchar* _map = map + mapstep*i + 1;
        _map[-1] = _map[src.cols] = 1;

        int* _mag = mag_buf[1] + 1; // take the central row
        ptrdiff_t magstep1 = mag_buf[2] - mag_buf[1];
        ptrdiff_t magstep2 = mag_buf[0] - mag_buf[1];

        const short* _x = dx.ptr<short>(i-1);
        const short* _y = dy.ptr<short>(i-1);

        if ((stack_top - stack_bottom) + src.cols > maxsize)
        {
            int sz = (int)(stack_top - stack_bottom);
            maxsize = std::max(maxsize * 3/2, sz + src.cols);
            stack.resize(maxsize);
            stack_bottom = &stack[0];
            stack_top = stack_bottom + sz;
        }

        int prev_flag = 0;
        for (int j = 0; j < src.cols; j++)
        {
#define CANNY_SHIFT 15
            const int TG22 = (int)(0.4142135623730950488016887242097*(1<<CANNY_SHIFT) + 0.5);

            int m = _mag[j];

            if (m > low)
            {
                int xs = _x[j];
                int ys = _y[j];
                int x = std::abs(xs);
                int y = std::abs(ys) << CANNY_SHIFT;

                int tg22x = x * TG22;

                if (y < tg22x)
                {
                    if (m > _mag[j-1] && m >= _mag[j+1]) goto __ocv_canny_push;
                }
                else
                {
                    int tg67x = tg22x + (x << (CANNY_SHIFT+1));
                    if (y > tg67x)
                    {
                        if (m > _mag[j+magstep2] && m >= _mag[j+magstep1]) goto __ocv_canny_push;
                    }
                    else
                    {
                        int s = (xs ^ ys) < 0 ? -1 : 1;
                        if (m > _mag[j+magstep2-s] && m > _mag[j+magstep1+s]) goto __ocv_canny_push;
                    }
                }
            }
            prev_flag = 0;
            _map[j] = uchar(1);
            continue;
__ocv_canny_push:
            if (!prev_flag && m > high && _map[j-mapstep] != 2)
            {
                CANNY_PUSH(_map + j);
                prev_flag = 1;
            }
            else
                _map[j] = 0;
        }

        // scroll the ring buffer
        _mag = mag_buf[0];
        mag_buf[0] = mag_buf[1];
        mag_buf[1] = mag_buf[2];
        mag_buf[2] = _mag;
    }

    // now track the edges (hysteresis thresholding)
    while (stack_top > stack_bottom)
    {
        uchar* m;
        if ((stack_top - stack_bottom) + 8 > maxsize)
        {
            int sz = (int)(stack_top - stack_bottom);
            maxsize = maxsize * 3/2;
            stack.resize(maxsize);
            stack_bottom = &stack[0];
            stack_top = stack_bottom + sz;
        }

        CANNY_POP(m);

        if (!m[-1])         CANNY_PUSH(m - 1);
        if (!m[1])          CANNY_PUSH(m + 1);
        if (!m[-mapstep-1]) CANNY_PUSH(m - mapstep - 1);
        if (!m[-mapstep])   CANNY_PUSH(m - mapstep);
        if (!m[-mapstep+1]) CANNY_PUSH(m - mapstep + 1);
        if (!m[mapstep-1])  CANNY_PUSH(m + mapstep - 1);
        if (!m[mapstep])    CANNY_PUSH(m + mapstep);
        if (!m[mapstep+1])  CANNY_PUSH(m + mapstep + 1);
    }

    // the final pass, form the final image
    const uchar* pmap = map + mapstep + 1;
    uchar* pdst = dst.ptr();
    for (int i = 0; i < src.rows; i++, pmap += mapstep, pdst += dst.step)
    {
        for (int j = 0; j < src.cols; j++)
            pdst[j] = (uchar)-(pmap[j] >> 1);
    }

    return dst;
}
