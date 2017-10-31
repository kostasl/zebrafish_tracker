#ifndef FISHMODEL_H
#define FISHMODEL_H

#include <string>

#include <cvBlob/cvblob.h>
#include <QDebug>
//#include <QApplication>


#include "ellipse_detect.h"
#include "zfttracks.h"


//#include "larvatrack.h" /If included here it causes circular search if fishModel Defs.

extern float gDisplacementThreshold;
extern int gFishTailSpineSegmentLength;
extern int gFishTailSpineSegmentCount;

/// \brief defines points along our custom linear spline that is fitted along the fish contour
typedef struct
{
    float x;
    float y; ///Position of Joint In Global Coordinates
    float angleRad;/// In Rads
} splineKnotf;

typedef std::vector<splineKnotf> t_fishspline;


typedef cv::KeyPoint zftblob;
typedef std::vector<zftblob> zftblobs;


const cv::Scalar TRACKER_COLOURMAP[] ={
                                CV_RGB(250,250,0), //Yellow
                                CV_RGB(255,0,0),
                                 CV_RGB(200,100,100),
                                 CV_RGB(150,200,50),
                                 CV_RGB(50,250,00),
                                 CV_RGB(150,150,00),
                                 CV_RGB(250,250,00),
                                 CV_RGB(200,200,80),
                                 CV_RGB(20,200,180),
                                 CV_RGB(200,20,180)};



/// \var typedef std::list<CvPoint2D64f> CvTrackPoints
/// \brief stores the stacked List of past centroid points that define this track.
/// \see CvPoint2D64f
/// \see CvTrack
typedef std::vector<cv::Point2f> CvTrackPoints;


class fishModel
{


public:
  fishModel();
  ~fishModel();
  //fishModel(cvb::CvTrack* track,cvb::CvBlob* blob);
  fishModel(zftblob blob,int bestTemplateOrientation,cv::Point ptTemplateCenter);

  void updateState(zftblob* fblob,double templatematchScore,int Angle, cv::Point2f bcentre,unsigned int nFrame);
  ///\note The lowest point in a rectangle is 0th vertex, and 1st, 2nd, 3rd vertices follow clockwise.
  /// Height is distance between 0th & 1st  (or 2nd & 3rd) vertices. And width is distance between 1st  & 2nd (or 0th & 3rd) vertices.
  /// Angle is calculated from the horizontal to the first edge of rectangle, counter clockwise.
  ///  Angle varies between -0 to -90 (I am unsure, what is the decisive factor of -0 or -90)
  float leftEyeAngle();
  float rightEyeAngle();
  float vergenceAngle();
  void resetSpine();
  void calcSpline(t_fishspline& outspline);
  void getSplineParams(t_fishspline& inspline,std::vector<double>& outparams);
  void setSplineParams(t_fishspline& inspline,std::vector<double>& inparams);

  cv::Point2f getPointAlongSpline(float z,t_fishspline& pspline);

  double distancePointToSpline(cv::Point2f ptsrc,t_fishspline& outspline);
  double getdeltaSpline(t_fishspline inspline, t_fishspline& outspline,int idxparam,double sgn);///
  //double fitSpineToContour(std::vector<std::vector<cv::Point> >& contours_body,int idxInnerContour,int idxOuterContour);
  double fitSpineToContour(cv::Mat& frameImg_grey, std::vector<std::vector<cv::Point> >& contours_body,int idxInnerContour,int idxOuterContour);
  void fitSpineToIntensity(cv::Mat &imgframeIn); //Uses Image Intensity Local Max to fit spline
  void GioGet_tailSpine(cv::Mat &src, cv::Point2i start, cv::Point2d tgt_start, int step_size, std::vector<cv::Point2i>& anchor_pts);


  void drawSpine(cv::Mat& outFrame);
  void drawBodyTemplateBounds(cv::Mat& outframe);
  friend std::ostream& operator<<(std::ostream& out, const fishModel& h);
  friend QTextStream& operator<<(QTextStream& out, const fishModel& h);

  zftID ID; /// Uid Of this Fish Instance

  cvb::CvLabel blobLabel; //Legacy BlobLabel

  std::vector<cv::Point> contour;
  std::vector<cv::Point> coreHull; /// core Body shape- no tail
  std::vector<cv::Point> coreTriangle; /// Core Body triangle Approximation


  unsigned int nCurrentFrame; ///<-Holds the frame Number of the last State Update
  double templateScore; //FishLIke Score - How well the detected model fish looks/matches the convolution of a fish template
  double leftEyeTheta; /// Theta is In Degrees
  double rightEyeTheta;/// Theta is In Degrees
  double bearingRads; /// Rads
  float bearingAngle; /// Theta is In Degrees


  tDetectedEllipsoid    leftEye;
  tDetectedEllipsoid    rightEye;
  cv::RotatedRect       bodyRotBound;

  cv::Point mouthPoint;
  cv::Point midEyePoint;
  cv::Point ptRotCentre; //Template Matching Body Centre
  cv::Point tailTopPoint;

  //cvb::CvTrack* track; ///Pointer to Track Structure containing motion - Note track has the same Id as this Fish
  //CvTrackPoints trackPointStack; /// <Holds list of past centroid positions along the track

  zftTrack zTrack;

  zftblob  zfishBlob; //Copy To assigned Blob structure
  t_fishspline spline; ///X-Y Coordinates of Fitted spline to contour

  int c_spineSegL;
  const int c_spinePoints   = gFishTailSpineSegmentCount;
  const int c_spineParamCnt = c_spinePoints+2;
private:

  //std::vect mmor<double> splineTheta; ///Angles of fitted Spine Points
};


///
/// \brief operator << //Overloaded Stream Operator // Output Current State Of The Track
/// \param out
/// \param h
/// \return
///
std::ostream& operator<<(std::ostream& out, const zftTrack& h);
QTextStream& operator<<(QTextStream& out, const zftTrack& h);
std::ostream& operator<<(std::ostream& out, const fishModel& h);
QTextStream& operator<<(QTextStream& out, const fishModel& h);


///Auxiliary


///
/// \brief fishModels list of model structures describing each visible fish
/// this list is maintained along with tracks - ie deletion/creation is done via matching to
/// blobs
///
typedef std::map<cvb::CvLabel,fishModel* > fishModels;


class CompareFishScore {
    public:
    bool operator()(fishModel*& t1, fishModel*& t2) // Returns true if t1 is greater than t2 /Ordering Highest 1st
    {
       return t1->templateScore > t2->templateScore;
    }
};

///
/// \brief qfishModels used to rank candidate fishModels according to match score
/// so as to obtain the best match among all candidate fishModels found, and remove the rest.
///
typedef std::priority_queue<fishModel*,std::vector<fishModel*>,CompareFishScore> qfishModels;



// /// \var typedef std::pair<CvID, fishModel *> IDFishModel pair for insertion into map list of fish
// /// /// \brief Pair (identification number, fishModel).
// /// \see CvID
// /// \see CvTrack
typedef std::pair<zftID, fishModel* > IDFishModel;




/// \fn inline void cvReleaseFishModes(fishModels &fishes)
/// \brief Clear Fish LIst
/// \param fishmodles List
/// \see
inline void ReleaseFishModels(fishModels &fishes)
{
  for (fishModels::iterator it=fishes.begin(); it!=fishes.end(); ++it)
  {
      fishModel* fish = (*it).second;
        //Let ReleaseTracks Handle This
//      if (fish->track)
//      {
//         fish->track->pointStack.clear();
//         delete fish->track;
//      }

      delete fish;
  }
  fishes.clear();
}



#endif // FISHMODEL_H
