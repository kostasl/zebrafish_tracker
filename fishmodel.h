#ifndef FISHMODEL_H
#define FISHMODEL_H

#include <string>
#include <cvBlob/cvblob.h>


class fishModel
{


public:
  fishModel();
  fishModel(cvb::CvTrack* track);

  float leftEyeAngle();
  float rightEyeAngle();
  float vergenceAngle();
  void resetSpine();
  void getSpline(std::vector<cv::Point2f>& outspline);
  double fitSpineToContour(std::vector<std::vector<cv::Point> >& contours_body,int idxInnerContour,int idxOuterContour);

  cvb::CvLabel blobLabel;
  cvb::CvID ID; /// Same as the track ID
  std::vector<cv::Point> contour;

  ///\note The lowest point in a rectangle is 0th vertex, and 1st, 2nd, 3rd vertices follow clockwise.
  /// Height is distance between 0th & 1st  (or 2nd & 3rd) vertices. And width is distance between 1st  & 2nd (or 0th & 3rd) vertices.
  /// Angle is calculated from the horizontal to the first edge of rectangle, counter clockwise.
  ///  Angle varies between -0 to -90 (I am unsure, what is the decisive factor of -0 or -90)
  cv::RotatedRect leftEyeRect;
  cv::RotatedRect rightEyeRect;
  std::vector<cv::Point> rightEyeHull;
  std::vector<cv::Point> leftEyeHull;
  std::vector<cv::Point> coreHull; /// core Body shape- no tail
  std::vector<cv::Point> coreTriangle; /// Core Body triangle Approximation

  double leftEyeTheta;
  double rightEyeTheta;
  double bearingRads;

  cv::Point leftEyePoint; /// Rotation Angle against Fish's head Midline
  cv::Point rightEyePoint;
  cv::Point mouthPoint;
  cv::Point midEyePoint;
  cv::Point tailTopPoint;

  cvb::CvTrack* track; ///Pointer to Track Structure containing motion - Note track has the same Id as this Fish

private:
  static int c_spinePoints  = 8;
  static int c_spineSegL    = 5;
  std::vector<cv::Point2f> spline; ///X-Y Coordinates of Fitted spline to contour
  std::vector<double> splineTheta; ///Angles of fitted Spine Points
};

#endif // FISHMODEL_H
