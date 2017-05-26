#ifndef LTROI_H
#define LTROI_H

//#include <cv.h> //Deprecated Header
#include <opencv2/core/core.hpp>

/// \typedef ltROI Custom Circular ROIs
/// \brief  Tracker Region of interest - used to define vials
typedef struct Circle
{
    cv::Point centre;
    unsigned int radius;

public:
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


    Circle(cv::Point a, cv::Point b)
    {
       centre = a;
       radius = cv::norm(b - centre);
    }

    bool contains(cv::Point pt)
    {
        return (cv::norm(pt - centre) < radius);
    }

    bool operator ==(const Circle& c2)
    {
        return (centre == c2.centre && radius == c2.radius);
    }
} ltROI; //Object to define region of interest

/// \typedef ltROIlist
///
typedef std::vector<ltROI> ltROIlist;

void addROI(ltROI& newRoi);
void deleteROI(cv::Point mousePos);
/// \fn Draw the Regions on the frame image
void drawROI();

/// \fn  ltROI* ltGetFirstROIContainingPoint(ltROIlist& vRoi, cv::Point pnt)
/// \brief Loop Over ROI
/// \params vRoi list
/// \params pnt to check membership for
ltROI* ltGetFirstROIContainingPoint(ltROIlist& vRoi ,cv::Point pnt);


#endif // LTROI_H
