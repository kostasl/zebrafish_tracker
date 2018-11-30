#ifndef ZFTTRACKS
#define ZFTTRACKS


#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d.hpp>
#include <string>
#include <QTextStream>

#define CV_TRACK_RENDER_ID            0x0001 ///< Print the ID of each track in the image. \see cvRenderTrack
#define CV_TRACK_RENDER_BOUNDING_BOX  0x0002 ///< Draw bounding box of each track in the image. \see cvRenderTrack
#define CV_TRACK_RENDER_TO_LOG        0x0010 ///< Print track info to log out. \see cvRenderTrack
#define CV_TRACK_RENDER_TO_STD        0x0020 ///< Print track info to log out. \see cvRenderTrack
#define CV_TRACK_RENDER_PATH          0x0100 ///< Draw polyline of track positions \see cvRenderTrack
#define CV_TRACK_RENDER_HIGHLIGHT     0x1000 ///Draw Colour Indicate Highlight

/// \brief Type of identification numbers.
typedef unsigned int zftID;
typedef unsigned int zfdID;


/// \var typedef std::vector<cv:Point2f> TrackPoints
/// \brief stores the stacked List of past centroid points that define this track.
/// \see CvPoint2D64f
/// \see CvTrack
typedef cv::Point2f zftTrackPoint;
typedef std::vector<zftTrackPoint> zftTrackPoints;
typedef std::vector<cv::Point> zftRenderedTrackPoints; //Used only For Rendering on image /Shorter vector And Integer based


///// \var typedef std::map<CvID, CvTrack *> CvTracks
///// \brief List of tracks.
///// \see CvID
///// \see CvTrack
//typedef std::map<CvID, CvTrack *> CvTracks;

///// \var typedef std::pair<CvID, CvTrack *> CvIDTrack
///// \brief Pair (identification number, track).
///// \see CvID
///// \see CvTrack
//typedef std::pair<CvID, CvTrack *> CvIDTrack;

//CvFont* defaultFont = NULL;



/// \brief Struct that contain information about one track.
/// \see CvID
/// \see ltROI
/// \see pointStack
/// \note zftTrack Maintains two list of track Points - the detailed floating point one is used for
/// analysis of movement, while the RenderPoints is optimized for rendering tracks, its integer and only holds
/// points which are displaced above D>=0.5, making the list shorted and quicker to render.
struct zftTrack
{
    zftTrack()
    {
        //Random colour
        int c1 =  rand() % 200 + 30;
        int c2 =  rand() % 200 + 30;
        int c3 =  rand() % 200 + 30;
        colour      = CV_RGB(c1,c2,c3);
        lifetime    = 0;
        active      = 0;
        inactive    = 0;
        effectiveDisplacement = 0.0;

    }

    //Default Constructor
    zftTrack(zftID ID):zftTrack()
    {
        id = ID;
    }

  zftID id; ///< Track identification number.
  //ltROI* pROI; ///< Pointer To Region of Interest structure to which this track belongs
  cv::Scalar colour = CV_RGB(255., 0., 0.); ///> Colourwhen drawing Countour

  cv::Rect boundingBox;
  cv::Point2f centroid; ///< Centroid.
  zftTrackPoints pointStack; /// <Holds list of past centroid positions along the track
  zftRenderedTrackPoints pointStackRender; //List Of Int Points Used for rendering Only

  double effectiveDisplacement; ///< Used to indicate a px speed measure so as to estimate possible blob distance from track on next frame.
  unsigned int lifetime; ///< Indicates how much frames the object has been in scene.
  unsigned int active; ///< Indicates number of frames that has been active from last inactive period.
  unsigned int inactive; ///< Indicates number of frames that has been missing.



};



/// \brief Render A zftracker Track
///
void zftRenderTrack(zftTrack& track, const cv::Mat& frameIn, cv::Mat& frameOut, unsigned short mode, int fontface,float fontScale );



std::ostream& operator<<(std::ostream& out, const zftTrack& h);
QTextStream& operator<<(QTextStream& out, const zftTrack& h);
#endif // ZFTTRACKS

