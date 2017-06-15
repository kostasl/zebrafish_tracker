#ifndef LARVATRACK_H
#define LARVATRACK_H


#include <iostream>
#include <sstream>
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

#include <GUI/mainwindow.h>

/// \file larvatrack.h
/// \brief OpenCV based zebrafish tracker header file.


class MainWindow;



struct fishModel
{

    cvb::CvLabel blobLabel;
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

    fishModel()
    {
        coreTriangle.push_back(cv::Point());
        coreTriangle.push_back(cv::Point());
        coreTriangle.push_back(cv::Point());
    }

};



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
/// \param maskedGrayImg
void detectZfishFeatures(cv::Mat& maskedGrayImg);
void checkPauseRun(MainWindow& win,int& keyboard,unsigned int nFrame);
bool saveImage(std::string frameNumberString,QString dirToSave,cv::Mat& img);
int countObjectsviaContours(cv::Mat& srcimg );
int processBlobs(cv::Mat& srcimg,cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString outFileCSV,std::string& frameNumberString,double& dMeanBlobArea);

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
int findContourClosestToPoint(std::vector<std::vector<cv::Point> >& contours,
                              std::vector<cv::Vec4i>& hierarchy,
                              cv::Point pt,
                              int level,
                              std::vector<std::vector<cv::Point> >& outhulls,
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
/// This is to recover Background substraction errors
/// \param frameImg
/// \param maskFGImg
///
void enhanceFishMask(cv::Mat& frameImg, cv::Mat& maskFGImg);

#endif // LARVATRACK_H
