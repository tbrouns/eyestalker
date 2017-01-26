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

#include "constants.h"
#include "parameters.h"
#include "structures.h"

void drawHaarDetector(cv::Mat&, int, int, int, int, cv::Vec3b);
void drawEdges(cv::Mat&, const std::vector<int>&, int, int, int, const cv::Vec3b&);
void drawOutline(cv::Mat&, const std::vector<std::vector<int> > &, int, int, int, const cv::Vec3b&, const cv::Vec3b&, const cv::Vec3b&);
void drawEllipse(cv::Mat&, const std::vector<double>&, int, int, int, int, const cv::Vec3b&);
void drawEllipseCross(cv::Mat&, double, double, const cv::Vec3b&);
void drawAll(cv::Mat&, detectionProperties);
void drawROI(cv::Mat&, int, int, int, int, cv::Vec3b);
