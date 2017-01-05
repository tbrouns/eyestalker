#ifndef CANNY
#define CANNY

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

cv::Mat cannyEdgeDetection(cv::Mat src, double low_thresh, double high_thresh, int aperture_size, bool L2gradient);

#endif // CANNY

