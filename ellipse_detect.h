#ifndef ELLIPSE_DETECT
#define ELLIPSE_DETECT

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
//#include <opencv2/bgsegm.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/video/background_segm.hpp>


typedef struct tDetectedEllipsoid{
    tDetectedEllipsoid(cv::RotatedRect r):rectEllipse(r){}

    cv::RotatedRect rectEllipse;
    int firscore;
} tDetectedEllipsoid;

typedef std::vector<tDetectedEllipsoid> tEllipsoids;

typedef struct tEllipsoidEdge {
    tEllipsoidEdge(cv::Point2f pt):ptEdge(pt) {}

    cv::Point2f ptEdge;
    int minorAxisLength;
} tEllipsoidEdge;

typedef std::vector<tEllipsoidEdge> tEllipsoidEdges;

int detectEllipses(cv::Mat& imgIn,cv::Mat& imgOut,tEllipsoids& vellipses);
void getEdgePoints(cv::Mat& imgEdgeIn,tEllipsoidEdges& vedgepoint);

#endif // ELLIPSE_DETECT

