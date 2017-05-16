#ifndef LARVATRACK_H
#define LARVATRACK_H


#include <iostream>
#include <sstream>



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
unsigned int processVideo(MainWindow& window_main, QString videoFilename,QString outFileCSV,unsigned int istartFrame);
void checkPauseRun(MainWindow& win,int& keyboard,std::string frameNumberString);
bool saveImage(std::string frameNumberString,QString dirToSave,cv::Mat& img);
int countObjectsviaContours(cv::Mat& srcimg );
int countObjectsviaBlobs(cv::Mat& srcimg,cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString outFileCSV,std::string& frameNumberString,double& dMeanBlobArea);

int saveTracks(cvb::CvTracks& tracks,QString filename,std::string frameNumber);
int saveTrackedBlobs(cvb::CvBlobs& blobs,QString filename,std::string frameNumber,cv::Rect& roi);
int saveTrackedBlobsTotals(cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString filename,std::string frameNumber,cv::Rect& roi);


void CallBackFunc(int event, int x, int y, int flags, void* userdata); //Mouse Callback


#endif // LARVATRACK_H
