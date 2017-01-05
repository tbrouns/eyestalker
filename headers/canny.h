#ifndef CANNY
#define CANNY

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

void cannyEdgeDetection(cv::Mat src, cv::Mat dst, double low_thresh, double high_thresh, int aperture_size, bool L2gradient);

#endif // CANNY

