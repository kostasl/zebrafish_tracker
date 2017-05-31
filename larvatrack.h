#ifndef LARVATRACK_H
#define LARVATRACK_H


#include <iostream>
#include <sstream>
#include <string>
#include <QDebug>

#include <QString>
#include <QApplication>
#include <QQmlApplicationEngine>
#include <QDir>
#include <QFileDialog>
#include <QTextStream>
#include <QElapsedTimer>
#include <QPixmap>

#include <cvBlob/cvblob.h>
#include <ltROI.h> //Defines the ROI types

#include <GUI/mainwindow.h>

class MainWindow;

/// \file larvatrack.h
/// \brief OpenCV based in-vial larva tracker header file.



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
int countObjectsviaBlobs(cv::Mat& srcimg,cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString outFileCSV,std::string& frameNumberString,double& dMeanBlobArea);

int saveTracks(cvb::CvTracks& tracks,QString filename,std::string frameNumber);
int saveTrackedBlobs(cvb::CvBlobs& blobs,QString filename,std::string frameNumber,ltROI& roi);
int saveTrackedBlobsTotals(cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString filename,std::string frameNumber,ltROI& roi);


void CallBackFunc(int event, int x, int y, int flags, void* userdata); //Mouse Callback


#endif // LARVATRACK_H
