#ifndef LARVATRACK_H
#define LARVATRACK_H


#include <config.h>
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

//#include <cvBlob/cvblob.h>
#include <ltROI.h> //Defines the ROI types
#include <fishmodel.h>
#include <foodmodel.h>
#include <zfttracks.h>
#include <GUI/mainwindow.h>

/// For Mem Usage //
#include <unistd.h>
#include <ios>
#include <fstream>


#undef _DEBUG

// Tail fitting, Use a complementary method to intensity to fit a spine to contour
// Useful as it can detect large deviations of tail from fish contour
#undef _USEFITSPINETOCONTOUR
//Run Tail Fit at as some interval, so as to check on the correctness of fit and reset if we have to
#define _USEPERIODICSPINETOCONTOUR_TEST //Periodic Scoring Of Tail



/// \file larvatrack.h
/// \brief OpenCV based zebrafish tracker header file.


class MainWindow;
class fishModel;
class vfish;

typedef cv::KeyPoint zftblob;
typedef std::vector<zftblob> zftblobs;

//typedef cv::KeyPoint zftblob;
//typedef std::vector<zftblob> zftblobs;



/// \fn processVideo
/// \brief  Main Loop function runs for each video in a list assumes 1st frame id is istartFrame.
/// \returns Returns the last frame number
/// \param videoFilename
/// \param outFileCSV
/// \param istartFrame
unsigned int processVideo(cv::Mat& fgMask,MainWindow& window_main, QString videoFilename,QFile& outdatafile,unsigned int istartFrame,unsigned int istopFrame);
//unsigned int getBGModelFromVideo(cv::Mat& fgMask,MainWindow& window_main,QString videoFilename,QString outFileCSV,unsigned int startFrameCount);
unsigned int trackImageSequencefiles(MainWindow& window_main);
//unsigned int trackVideofiles(MainWindow& window_main,QString outputFile);
unsigned int trackVideofiles(MainWindow& window_main,QString outputFile,QStringList invideonames,unsigned int istartFrame,unsigned int istopFrame);

/// \fn processFrame - Process blob morphology, Extract features tracks
///
//void processFrame(cv::Mat& frame,cv::Mat& fgMask, unsigned int nFrame,cv::Mat& outframe,cv::Mat& outframeHead);
void processFrame(MainWindow& window_main,const cv::Mat& frame,cv::Mat& fgMask, unsigned int nFrame,cv::Mat& outframe,cv::Mat& outframeHeadEyeDetected,cv::Mat& frameHead);
void drawFrameText(MainWindow& window_main, uint nFrame,uint nLarva,uint nFood,cv::Mat& outframe);
//bool updateBGFrame(cv::Mat& frame,cv::Mat& fgMask, unsigned int nFrame);

///
/// \brief detectZfishFeatures - Used to create geometric representations of main zebrafish Features : Eyes, Body, tail
/// these are saved as point arrays on which angles and other measurements can be obtained
/// \param fullImg - Raw Img captured
/// \param fullImgOut - The labels and output Produced by the Fnct is written on this canvas
/// \param  headImgOut Just the head INset
/// \param Mask with Fish oNly FG
/// \param main inner and extrernal fish Contours for each fish
/// \param 1 level hierarchy of contours (outer inner)
void detectZfishFeatures(MainWindow& window_main,const cv::Mat& fullImgIn,cv::Mat& fullImgOut,cv::Mat& headImgOut,cv::Mat& outheadImgOutProcessed, cv::Mat& maskfishFGImg, std::vector<std::vector<cv::Point> >& contours_body,std::vector<cv::Vec4i>& hierarchy_body);
///
/// \brief UpdateFishModels Use Tracks  to update persistent fishModels
/// \param vfishmodels
/// \param fishtracks
//void UpdateFishModels(cv::Mat& fullImgIn,fishModels& vfishmodels,cvb::CvTracks& fishtracks); /// Deprecated
///Updated Version With Blobs
void UpdateFishModels(const cv::Mat& fullImgIn,fishModels& vfishmodels,zftblobs& fishtracks,unsigned int nFrame,cv::Mat& frameOut);
void UpdateFoodModels(const cv::Mat& maskedImg_gray,foodModels& vfoodmodels,zfdblobs& foodblobs,unsigned int nFrame,bool bAllowNew); //,cv::Mat& frameOut
int processFoodOpticFlow(const cv::Mat frame_grey,const cv::Mat frame_grey_prev,foodModels& vfoodmodels,unsigned int nFrame,zftblobs& vPreyKeypoints_next  );

bool checkKeypointExistsAndRemove(zfdblobs& vfoodblobs_spare,zfdblob& FoodBlob); //Matches Keypoints/Blobs based on position and removes match from vector

void checkPauseRun(MainWindow* win,int keyboard,unsigned int nFrame);
void keyCommandFlag(MainWindow* win, int keyboard,unsigned int nFrame);
bool saveImage(QString frameNumberString,QString dirToSave,QString filenameVid,cv::Mat& img);
int countObjectsviaContours(cv::Mat& srcimg );
//int processBlobs(IplImage* srcframeImg,cv::Mat& maskimg,cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString outDirCSV,std::string& frameNumberString,double& dMeanBlobArea);

/// Updated Blobs Detector- Fish Specific
int processFishBlobs(cv::Mat& frame,cv::Mat& maskimg,cv::Mat& frameOut,zftblobs& ptFishblobs);
int processFoodBlobs(const cv::Mat& frame,const cv::Mat& maskimg,cv::Mat& frameOut,zftblobs& ptFoodblobs);

int saveTracks(fishModels& vfish,foodModels& vfood, QFile& data, QString frameNumber);
int saveFoodTracks(fishModels& vfish,foodModels& vfood, QFile& fooddata,QString frameNumber);
bool openDataFile(QString filepathCSV,QString filenameVid,QFile& data,QString strpostfix="_tracks");
void closeDataFile(QFile& data);
void removeDataFile(QFile& data);
bool resetDataRecording(QFile& outdatafile,QString strpostfix); //Uses Global File Info

void writeFishDataCSVHeader(QFile& data);
void writeFoodDataCSVHeader(QFile& data);

//int saveTrackedBlobs(cvb::CvBlobs& blobs,QString filename,std::string frameNumber,ltROI& roi);
//int saveTrackedBlobsTotals(cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString filename,std::string frameNumber,ltROI& roi);

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

///Simpler One Using Only distance from Point
int findMatchingContour(std::vector<std::vector<cv::Point> >& contours,
                              std::vector<cv::Vec4i>& hierarchy,
                              cv::Point pt,
                              int level);

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
//void enhanceFishMask(cv::Mat& frameImg, cv::Mat& maskFGImg,std::vector<std::vector<cv::Point> >& fishbodycontours ,std::vector<cv::Vec4i>& fishbodyhierarchy);
void enhanceMask(const cv::Mat& frameImg, cv::Mat& maskFGImg,cv::Mat& outFishMask,cv::Mat& outFoodMask,std::vector<std::vector<cv::Point> >& fishbodycontours, std::vector<cv::Vec4i>& fishbodyhierarchy);

///
/// \brief findIndexClosesttoPoint - Naive Nearest neighbour finder
/// \param vPointChain array of points
/// \param pt reference point
/// \return index in array of closest point
///
int findIndexClosesttoPoint(std::vector<cv::Point> vPointChain,cv::Point pt);

/// \brief locates tail by taking curve derivaties
///
int smoothContour(const cv::Mat& frameImg, cv::Mat& fgMask,std::vector<cv::Point>& curve);

/// \brief Returns the index of the point furthest away from the provided tail point idx on the curve
/// That furthest from the tail position is likely the larva's head
int findAntipodePointinContour(int idxTail, std::vector<cv::Point>& curve,cv::Point ptCentroid, cv::Point& ptHead,cv::Point& ptTail);


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

///
/// \brief process_mem_usage Attempts a read of MemUsage
/// \param vm_usage
/// \param resident_set
///
void process_mem_usage(double& vm_usage, double& resident_set);

//////  Stream Output Operators ////


/// Auxiliriary Functions
template < typename T > std::string to_string( const T& n )
   {
       std::ostringstream stm ;
       stm << n ;
       return stm.str() ;
   }




///
/// \brief opencv event Callback functions
///
void CallBackFunc(int event, int x, int y, int flags, void* userdata); //Mouse Callback
void CallBackHistFunc(int event, int x, int y, int flags, void* userdata); //Mouse Callback For the Histogram Window
void thresh_callback(int, void* );

/// To Empty Priority Queue
///
///
template <class T, class S, class C>
void clearpq2(std::priority_queue<T, S, C>& q){
    q=std::priority_queue<T, S, C>();
}



#endif // LARVATRACK_H

