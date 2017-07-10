#ifndef LARVATRACK_H
#define LARVATRACK_H


#include <iostream>
#include <sstream>
#include <iomanip> //for setprecision
#include <limits>
#include <string>
#include <QDebug>

#include <QString>
#include <QApplication>
//#include <QQmlApplicationEngine>
#include <QDir>
#include <QFileDialog>
#include <QTextStream>
#include <QElapsedTimer>
#include <QPixmap>

#include <cvBlob/cvblob.h>
#include <ltROI.h> //Defines the ROI types
#include <fishmodel.h>

#include <GUI/mainwindow.h>

/// \file larvatrack.h
/// \brief OpenCV based zebrafish tracker header file.


class MainWindow;



///
/// \brief fishModels list of model structures describing each visible fish
/// this list is maintained along with tracks - ie deletion/creation is done via matching to
/// blobs
///
typedef std::map<cvb::CvLabel,fishModel* > fishModels;


/// \var typedef std::pair<CvID, fishModel *> CvIDFishModel pair for insertion into map list of fish
/// /// \brief Pair (identification number, fishModel).
/// \see CvID
/// \see CvTrack
typedef std::pair<cvb::CvID, fishModel* > CvIDFishModel;



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



/// \fn processVideo
/// \brief  Main Loop function runs for each video in a list assumes 1st frame id is istartFrame.
/// \returns Returns the last frame number
/// \param videoFilename
/// \param outFileCSV
/// \param istartFrame
unsigned int processVideo(cv::Mat& fgMask,MainWindow& window_main, QString videoFilename,QString outFileCSV,unsigned int istartFrame);
unsigned int getBGModelFromVideo(cv::Mat& fgMask,MainWindow& window_main,QString videoFilename,QString outFileCSV,unsigned int startFrameCount);
unsigned int trackImageSequencefiles(MainWindow& window_main);
unsigned int trackVideofiles(MainWindow& window_main);
/// \fn processFrame - Process blob morphology, Extract features tracks
///
void processFrame(cv::Mat& frame,cv::Mat& fgMask,cv::Mat& frameMasked, unsigned int nFrame);
bool updateBGFrame(cv::Mat& frame,cv::Mat& fgMask, unsigned int nFrame);

///
/// \brief detectZfishFeatures - Used to create geometric representations of main zebrafish Features : Eyes, Body, tail
/// these are saved as point arrays on which angles and other measurements can be obtained
/// \param fullImg - Raw Img captured
/// \param Mask with Fish oNly FG
/// \param main inner and extrernal fish Contours for each fish
/// \param 1 level hierarchy of contours (outer inner)
void detectZfishFeatures(cv::Mat& fullImg,cv::Mat& maskfishFGImg,std::vector<std::vector<cv::Point> >& contours_body,std::vector<cv::Vec4i>& hierarchy_body);

///
/// \brief UpdateFishModels Use Tracks  to update persistent fishModels
/// \param vfishmodels
/// \param fishtracks
///
void UpdateFishModels(fishModels& vfishmodels,cvb::CvTracks& fishtracks);

void checkPauseRun(MainWindow* win,int keyboard,unsigned int nFrame);
void keyCommandFlag(MainWindow* win, int keyboard,unsigned int nFrame);
bool saveImage(std::string frameNumberString,QString dirToSave,cv::Mat& img);
int countObjectsviaContours(cv::Mat& srcimg );
int processBlobs(IplImage* srcframeImg,cv::Mat& maskimg,cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString outDirCSV,std::string& frameNumberString,double& dMeanBlobArea);

int saveTracks(cvb::CvTracks& tracks,QString filename,std::string frameNumber);
int saveTrackedBlobs(cvb::CvBlobs& blobs,QString filename,std::string frameNumber,ltROI& roi);
int saveTrackedBlobsTotals(cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString filename,std::string frameNumber,ltROI& roi);

void CallBackFunc(int event, int x, int y, int flags, void* userdata); //Mouse Callback
/**
* @function thresh_callback
*/
void thresh_callback(int, void* );

///
/// \brief findContourClosestToPoint Looks for the inner contour in a 2 level hierarchy that matches the point coords
/// \param contours
/// \param hierarchy
/// \param pt
/// \param level - The required hierarchy level description of the contour being searched for
/// \param OUT outhull
/// \param Out fittedEllipse - pointer to array of Rotated rect fitted ellipsoids
/// \return Index of *child*/Leaf contour closest to point
///
int findMatchingContour(std::vector<std::vector<cv::Point> >& contours,
                              std::vector<cv::Vec4i>& hierarchy,
                              cv::Point pt,
                              int level,
                              std::vector<cv::Point>& matchhull,
                              std::vector<cv::RotatedRect>& outfittedEllipse);



///
/// \brief fitfishCoreTriangle Sets a fixed position to represent fish features Guesses tail point
/// \param maskedfishFeature Image containing only the identified fish pixels
/// \param sfish
/// \param contours_body
/// \param idxInnerContour Pass Index for inner fish body contour (as segregated by morph on thresholded image
/// \param idxOuterContour Pass index of the outer whole fish contour
/// \return
///
bool fitfishCoreTriangle(cv::Mat& maskedfishFeature,fishModel& sfish,std::vector<std::vector<cv::Point> >& contours_body,int idxInnerContour,int idxOuterContour);


///
/// \brief enhanceFishMask Looks for fish countours and draws them onto the FG mask so as to enhance features
/// This is to recover Background substraction errors and obtain accurate Fish Masks before Blob Detection
/// \param frameImg
/// \param maskFGImg
///
void enhanceFishMask(cv::Mat& frameImg, cv::Mat& maskFGImg,std::vector<std::vector<cv::Point> >& fishbodycontours ,std::vector<cv::Vec4i>& fishbodyhierarchy);


///
/// \brief findIndexClosesttoPoint - Naive Nearest neighbour finder
/// \param vPointChain array of points
/// \param pt reference point
/// \return index in array of closest point
///
int findIndexClosesttoPoint(std::vector<cv::Point> vPointChain,cv::Point pt);

/// \brief Find point Furthest Along closed outline contour
/// Assume points belong to closed contour - find max distance between points
/// traversing clock and anticlockwise
int maxChainDistance(std::vector<cv::Point> vPointChain,int idx,int idy);



/// Find Eye Orientation by finding the angle of a fixed lentgh line that gives
///  maximum intensity
/// \param ptOutA, B two points defining the best  line fit
/// \return Best Ange In rads
double findEyeOrientation(cv::Mat& frameFish_gray, cv::Point2f& ptEyeCenter,std::vector<cv::Point>& outvEyeContour);


void makeEllipse(cv::Point2f ptcenter,double angle,double a, double b, std::vector<cv::Point>& voutEllipse);
/// Auxiliriary Functions

template < typename T > std::string to_string( const T& n )
   {
       std::ostringstream stm ;
       stm << n ;
       return stm.str() ;
   }

#endif // LARVATRACK_H

