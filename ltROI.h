#ifndef LTROI_H
#define LTROI_H

//#include <cv.h> //Deprecated Header
#include <opencv2/core/core.hpp>

/// \typedef ltROI
/// \brief Larva Tracker Region of interest - used to define vials
typedef cv::Rect ltROI; //Object to define region of interest
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
