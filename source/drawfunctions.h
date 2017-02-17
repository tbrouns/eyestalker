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

#include "constants.h"
#include "parameters.h"
#include "structures.h"

#include <iostream>

void drawEdges  (cv::Mat&, const AOIProperties&, const cv::Vec3b&, const std::vector<int>&);
void drawOutline(cv::Mat&, const AOIProperties&, const cv::Vec3b&, const std::vector<std::vector<int>>&);
void drawEllipse(cv::Mat&, const AOIProperties&, const cv::Vec3b&, const cv::Vec3b&, const cv::Vec3b&, const std::vector<double>&);
void drawAOI    (cv::Mat&, const AOIProperties&, const cv::Vec3b&);
void drawAll    (cv::Mat&, const drawVariables&);
void drawCross  (cv::Mat&, double, double, const cv::Vec3b&);
