#ifndef LTROI_H
#define LTROI_H

//#include <cv.h> //Deprecated Header
#include <vector>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/memory.hpp>


//Inherit and Add Serialization
class tROIpt: public cv::Point
{
    public:

    tROIpt():cv::Point(){};

    // This method lets cereal know which data members to serialize
      template<class Archive>
      void serialize(Archive & archive)
      {
        archive(CEREAL_NVP(x)); // serialize things by passing them to the archive
        archive(CEREAL_NVP(y)); // serialize things by passing them to the archive
      }

};



/// \typedef ltROI Custom  ROIs
/// \brief  Tracker Region of interest - used to define vials
typedef struct tRoi
{


    enum RoiType {Circle, Polygon};

    RoiType mType;
    std::vector<cv::Point> vPoints;

    cv::Point centre;

    cv::Scalar mColor = CV_RGB(255,5,5);
    unsigned int radius;

public:
    tRoi(); //Constructor

    inline int x()
    {
        return centre.x;
    }

    inline int y()
    {
        return centre.y;
    }

    inline unsigned int width()
    {
        return radius;
    }


    //Construct A Circular ROI given 2 points
    tRoi(cv::Point a, cv::Point b)
    {
       mType = RoiType::Circle;
       centre = a;
       vPoints.push_back(centre);
       radius = cv::norm(b - centre);
    }

    //Construct a polygon Roi Given A list of Points
    tRoi(std::vector<cv::Point>& pvPoints )
    {
        mType = RoiType::Polygon;
        vPoints = pvPoints; //local Copy of Points

        //Find the centre of this ROI
        cv::Moments mm = cv::moments(vPoints);

        centre.x = mm.m10/mm.m00;
        centre.y = mm.m01/mm.m00;


    }



   void  draw(cv::Mat& frame)
    {
        if (mType == RoiType::Circle)
        {
           cv::circle(frame,centre,radius,cv::Scalar(0,0,250),2);
        }

        if (mType == RoiType::Polygon)
        {
            //polylines(InputOutputArray img, InputArrayOfArrays pts, bool isClosed, const Scalar& color, int thickness=1, int lineType=8, int shift=0 )¶
            cv::polylines(frame,vPoints,true,mColor,1);
            //Draw Anchor points
            for (std::vector<cv::Point>::iterator it = vPoints.begin() ; it != vPoints.end(); ++it)
            {
                cv::circle(frame,*it,1,cv::Scalar(0,0,250),2);
            }

        }

    }


   void  drawMask(cv::Mat& frame)
    {
        if (mType == RoiType::Circle)
        {
           //cv::circle(frame,centre,radius,cv::Scalar(0,0,250),2);
           cv::circle(frame,centre,radius,CV_RGB(255,255,255),-1);
        }

        if (mType == RoiType::Polygon)
        {
            //void fillPoly(Mat& img, const Point** pts, const int* npts, int ncontours, const Scalar& color, int lineType=8, int shift=0, Point offset=Point() )¶
             //void fillConvexPoly(Mat& img, const Point* pts, int npts, const Scalar& color, int lineType=8, int shift=0)
            cv::fillConvexPoly(frame,vPoints,CV_RGB(255,255,255));
//            //Draw Anchor points
//            for (std::vector<cv::Point>::iterator it = vPoints.begin() ; it != vPoints.end(); ++it)
//            {
//                cv::circle(frame,*it,1,cv::Scalar(0,0,250),2);

//            }

        }

    }


 // Check If In Point Is Within ROI - Assuming Certain Size So as to Make A soft Boundary
   bool contains(cv::Point pt,double objSize = 1)
    {
        if (mType == RoiType::Circle)
            return (cv::norm(pt - centre) <= (radius+objSize));

        if (mType == RoiType::Polygon) //Check Distance On OUtside Is greater Than Obj Size - Only when objSize Param Is > 1
            return (cv::pointPolygonTest(vPoints,pt,(objSize > 1)) > -objSize );

        return false;
    }


    bool operator ==(const tRoi& c2)
    {
        if (c2.mType == RoiType::Circle)
            return (centre == c2.centre && radius == c2.radius);
    }

    // This method lets cereal know which data members to serialize
      template<class Archive>
      void serialize(Archive & archive)
      {
        archive(CEREAL_NVP(centre)); // serialize things by passing them to the archive
      }


} ltROI; //Object to define region of interest

//Need to register As archivable classes
//CEREAL_REGISTER_TYPE(ltROI);
//CEREAL_REGISTER_TYPE(tROIpt);


/// \typedef ltROIlist
///
typedef std::vector<ltROI> ltROIlist;

void addROI(ltROI& newRoi);
void deleteROI(cv::Point mousePos);
/// \fn Draw the Regions on the frame image
void drawAllROI(cv::Mat& frame);
void drawUserDefinedPoints(cv::Mat& frame);

/// \fn  ltROI* ltGetFirstROIContainingPoint(ltROIlist& vRoi, cv::Point pnt)
/// \brief Loop Over ROI
/// \params vRoi list
/// \params pnt to check membership for
ltROI* ltGetFirstROIContainingPoint(ltROIlist& vRoi ,cv::Point pnt,float objectRadius );
// Check if point is contained in any ROI found int the  global list of ROIs
bool pointIsInROI(cv::Point pt,float objectRadius ); //

#endif // LTROI_H
