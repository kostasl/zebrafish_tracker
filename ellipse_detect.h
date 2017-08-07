#ifndef ELLIPSE_DETECT
#define ELLIPSE_DETECT

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
//#include <opencv2/bgsegm.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/video/background_segm.hpp>

int detectEllipses(cv::Mat& imgIn,cv::Mat& imgOut,std::vector<cv::RotatedRect>& vellipses);
void getEdgePoints(cv::Mat& imgEdgeIn,std::vector<cv::Point>& vedgepoint);

#endif // ELLIPSE_DETECT

