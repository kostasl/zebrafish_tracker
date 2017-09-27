#ifndef ELLIPSE_DETECT
#define ELLIPSE_DETECT

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
//#include <opencv2/bgsegm.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/video/background_segm.hpp>



///\note Consistency between pts and RotRect is not checked
typedef struct tDetectedEllipsoid{
    //cv::RotatedRect(ptxy0,cv::Size2f(2*a,2*idx), alpha*(180/M_PI))
    tDetectedEllipsoid(cv::Point pt0,cv::Point pt1,cv::Point pt2,int score,cv::RotatedRect r):rectEllipse(r){


        ptAxisMj1 = pt1; //Major Axis Point 1;
        ptAxisMj2 = pt2; //Major Axis Point 2;
        fitscore = score;
    }

    tDetectedEllipsoid(cv::RotatedRect r,int score):rectEllipse(r){

        fitscore = score;
    }

    //Operator for Priority Ordering
//    bool operator<(const tDetectedEllipsoid& b) {
//      return this->fitscore < b.fitscore; //Max Heap
//    }

    cv::RotatedRect rectEllipse;
    int fitscore;
    cv::Point  ptAxisMj1;
    cv::Point  ptAxisMj2;

} tDetectedEllipsoid;


//Operator for Priority Ordering
bool operator<(const tDetectedEllipsoid& a,const tDetectedEllipsoid& b);

typedef std::vector<tDetectedEllipsoid> tEllipsoids;

typedef struct tEllipsoidEdge {
    tEllipsoidEdge(cv::Point2f pt):ptEdge(pt) {}

    cv::Point2f ptEdge;
    int minorAxisLength;
} tEllipsoidEdge;

typedef std::vector<tEllipsoidEdge> tEllipsoidEdges;

//int detectEllipses(cv::Mat& imgIn,cv::Mat& imgOut,tEllipsoids& vellipses);
//int detectEllipses(cv::Mat& imgIn,cv::Mat& imgOut,int angleDeg,tEllipsoids& vellipses);
int detectEllipse(tEllipsoidEdges& vedgePoints_all, std::priority_queue<tDetectedEllipsoid>& qEllipsoids);
///
/// \brief detectEllipses
/// \param pimgIn
/// \param imgEdge
/// \param imgOut
/// \param angleDeg
/// \param vellipses
/// \param outHeadFrameProc Return a close Up of the head with the detected shapes drawn on
/// \return
///
int detectEllipses(cv::Mat& pimgIn,cv::Mat imgEdge,cv::Mat& imgOut,int angleDeg,tEllipsoids& vellipses,cv::Mat& outHeadFrameProc);
void getEdgePoints(cv::Mat& imgEdgeIn,tEllipsoidEdges& vedgepoint);
void getEdgePoints(std::vector<cv::Point>& contour,tEllipsoidEdges& vedgepoint);

void show_histogram(std::string const& name, cv::Mat1b const& image);
#endif // ELLIPSE_DETECT

