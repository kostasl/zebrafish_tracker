#ifndef FISHMODEL_H
#define FISHMODEL_H

#include <string>
#include <cvBlob/cvblob.h>

class fishModel
{


public:
  fishModel();
  fishModel(cvb::CvTrack* track);

  cvb::CvLabel blobLabel;
  cvb::CvID ID; /// Same as the track ID
  std::vector<cv::Point> contour;
  cv::RotatedRect leftEyeHull;
  cv::RotatedRect rightEyeHull;
  std::vector<cv::Point> coreTriangle;

  cv::Point leftEyePoint; /// Rotation Angle against Fish's head Midline
  cv::Point rightEyePoint;
  double leftEyeTheta;
  double rightEyeTheta;
  cv::Point tailTopPoint;
  double bearingRads;
  cv::Point tailSplinePoints[8];

  cvb::CvTrack* track; ///Pointer to Track Structure containing motion - Note track has the same Id as this Fish


};

#endif // FISHMODEL_H
