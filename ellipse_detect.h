#ifndef ELLIPSE_DETECT
#define ELLIPSE_DETECT

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
//#include <opencv2/bgsegm.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/video/background_segm.hpp>
#include <queue>

//fwd declaration of vector of DetectedEllispods
//class tEllipsoids;

///\note Consistency between pts and RotRect is not checked
typedef struct tDetectedEllipsoid{
    tDetectedEllipsoid():rectEllipse()
    {
        fitscore = 0;
        cLabel = 'N'; //Not set / Default value
        ptAxisMj1.x = 0;
        ptAxisMj1.y = 0;
        ptAxisMj2.x = 0;
        ptAxisMj2.y = 0;
        nsamples = 1;

    }

    //Initialiazes object using mean values from vector of ellipsoids provided
    tDetectedEllipsoid(const std::vector<tDetectedEllipsoid>& vEll):tDetectedEllipsoid()
    {
        //Iterate all ellipsoids and obtain mean values

        //Note as position of MjAxis pt1 pt2 is not consistently lower/ upper
        // we need to avg in a directional way so the means are consistent for lower/upper point
        nsamples = vEll.size();
        for (int i=0;i<vEll.size();i++)
        {
            tDetectedEllipsoid rEll = vEll[i];
            if (rEll.ptAxisMj1.y > rEll.ptAxisMj2.y)
            { //ptAxisMj1 is lower than ptAxisMj2
                ptAxisMj1.y +=  rEll.ptAxisMj1.y ;
                ptAxisMj1.x +=  rEll.ptAxisMj1.x ;
                ptAxisMj2.y +=  rEll.ptAxisMj2.y ;
                ptAxisMj2.x +=  rEll.ptAxisMj2.x ;
            }
            else
           {
               ptAxisMj1.y +=  rEll.ptAxisMj2.y;
               ptAxisMj1.x +=  rEll.ptAxisMj2.x;
               ptAxisMj2.y +=  rEll.ptAxisMj1.y;
               ptAxisMj2.x +=  rEll.ptAxisMj1.x;
          }

            fitscore += rEll.fitscore;
        } //lOOP through all ellipsoids

        //Calc empirical mean /

        ptAxisMj1.y = ptAxisMj1.y / (float)nsamples;
        ptAxisMj1.x = ptAxisMj1.x / (float)nsamples;
        ptAxisMj2.y = ptAxisMj2.y / (float)nsamples;
        ptAxisMj2.x = ptAxisMj2.x / (float)nsamples;

        fitscore = fitscore/ (float)nsamples;

        cv::Point2f mjAxisLine = ptAxisMj2-ptAxisMj1;
        /// \todo set the other bounding rect points accordingly
        this->rectEllipse.angle = std::atan2(mjAxisLine.y,mjAxisLine.x) * 180.0/CV_PI+90.0;
    }

    //cv::RotatedRect(ptxy0,cv::Size2f(2*a,2*idx), alpha*(180/M_PI))
    /// \todo recalc all rect points / not just angle
    tDetectedEllipsoid(cv::Point2f pt0,cv::Point2f pt1,cv::Point2f pt2,int score,cv::RotatedRect r):rectEllipse(r)
    {
        if (pt1.y < pt2.y) //set so ptMj1 is the lower image point of the ellipse
        { //Pt1 is above pt2
            ptAxisMj1 = pt2; //Major Axis Point 1;
            ptAxisMj2 = pt1; //Major Axis Point 2;
        }
        else
        {//Pt1 is below pt2
            ptAxisMj1 = pt1; //Major Axis Point 1;
            ptAxisMj2 = pt2; //Major Axis Point 2;
        }

        fitscore = score;
        cv::Point2f mjAxisLine = ptAxisMj2-ptAxisMj1;
        r.angle = std::atan2(mjAxisLine.y,mjAxisLine.x) * 180.0/M_PI+90.0;
    }

    tDetectedEllipsoid(cv::RotatedRect r,int score):rectEllipse(r){

        if (r.angle >= 90)
            rectEllipse.angle = r.angle-180.0;
        if (r.angle <= -90)
            rectEllipse.angle = 180.0+r.angle;


        fitscore = score;
        ptAxisMj1.x = r.center.x + r.size.height*sin(-r.angle*M_PI/180.0)/3.0;
        ptAxisMj1.y = r.center.y + r.size.height*cos(-r.angle*M_PI/180.0)/3.0;

        ptAxisMj2.x = r.center.x - r.size.height*sin(-r.angle*M_PI/180.0)/3.0;
        ptAxisMj2.y = r.center.y - r.size.height*cos(-r.angle*M_PI/180.0)/3.0;

        //cv::Point2f ptBoxPts[4];
        //r.points(ptBoxPts);
        //ptAxisMj1 = r.center + (ptBoxPts[0]-(cv::Point2f)r.center)+(ptBoxPts[1]-ptBoxPts[0])/2.0;
        //ptAxisMj2 = r.center + (ptBoxPts[2]-(cv::Point2f)r.center); //+(ptBoxPts[2]-ptBoxPts[3])/2.0;
    }

    //Operator for Priority Ordering
//    bool operator<(const tDetectedEllipsoid& b) {
//      return this->fitscore < b.fitscore; //Max Heap
//    }

    cv::RotatedRect rectEllipse;
    int fitscore;
    int nsamples;
    cv::Point2f  ptAxisMj1;
    cv::Point2f  ptAxisMj2;
    float stdDev;
    char cLabel; ///A label char used to distinguish sets of detect ellipses (ie left right Eye here )



    //Returns corrected angle as used for reporting eye angles
    float getEyeAngle()
    {
        //float fEyeTheta;
        //fEyeTheta      = rectEllipse.angle;
//        if (fEyeTheta > 90)
//             fEyeTheta      = rectEllipse.angle-90;
//        if (fEyeTheta < -30)
//             fEyeTheta      = rectEllipse.angle+90;
    return (rectEllipse.angle);
    }


} tDetectedEllipsoid;

typedef std::vector<tDetectedEllipsoid> tEllipsoids;


//Operator for Priority Ordering
bool operator<(const tDetectedEllipsoid& a,const tDetectedEllipsoid& b);

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
/// \brief detectEllipses Search For Ellipsoids around the position of the eyes in the fish head isolated image
int detectEyeEllipses(cv::Mat& pimgIn,tEllipsoids& vLellipses,tEllipsoids& vRellipses,cv::Mat& outHeadFrameMonitor,cv::Mat& outHeadFrameProc);
// Uses Sampling around an arc below the eyes to determine Appropriate Eye Segmentation theshold (Max N values are used)
//int getEyeSegThreshold(cv::Mat& pimgIn,cv::Point2f ptcenter,std::vector<cv::Point>& ellipseSample_pts);
std::vector<int> getEyeSegThreshold(cv::Mat& pimgIn,cv::Point2f ptcenter,std::vector<cv::Point>& ellipseSample_pts,int& minVal,int& maxVal);
void getEdgePoints(cv::Mat& imgEdgeIn,tEllipsoidEdges& vedgepoint);
void getEdgePoints(std::vector<cv::Point>& contour,tEllipsoidEdges& vedgepoint);


/// Functions for optimizing ellipse detection to only consider connected pairs of points (ie points that belong to the same line/curve
void getEdgePoints(std::vector<cv::Point>& contour,tEllipsoidEdges& vedgepoint);/// Fills A list with  point coords where pixels (edges image) are above a threshold (non-zero)
void getPointsAlongEdge(cv::Mat imgEdgeIn,cv::Point2f startpt,tEllipsoidEdges& vedgepoint);
void getConnectedEdgePoints(cv::Mat& imgEdgeIn,cv::Point2f startpt,tEllipsoidEdges& vedgepoint);


void show_histogram(std::string const& name, cv::Mat1b const& image);
/// Draws LInes On Upsampled Head Image showing the major axis of ellipses
float drawExtendedMajorAxis(cv::Mat& outHeadFrameMonitor,tDetectedEllipsoid& ellEye,cv::Scalar col);
#endif // ELLIPSE_DETECT

