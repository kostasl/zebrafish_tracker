#ifndef FGMASKPROCESSING_H
#define FGMASKPROCESSING_H

#include <opencv2/opencv.hpp>
#include <GUI/mainwindow.h>
#include <QString>
#include <QApplication>

unsigned int getBGModelFromVideo(cv::Mat& fgMask,MainWindow& window_main,QString videoFilename,QString outFileCSV,unsigned int MOGhistoryLength);
bool updateBGFrame(cv::Mat& frame, cv::Mat& fgMask, unsigned int nFrame,uint MOGhistory);
void extractFGMask(cv::Mat& frameImg_gray,cv::Mat fgStaticMaskIn,cv::Mat& fgMaskInOut,cv::Mat& fgFrameOut,double dLearningRate);
int getMaxInflectionAndSmoothedContour(std::vector<cv::Point>& curve); // Curve Processing to find Tail/maxInflection point and simplify contour
void enhanceMasks(const cv::Mat& frameImg, cv::Mat& fgMask,cv::Mat& outFishMask,cv::Mat& outFoodMask,std::vector<std::vector<cv::Point> >& outfishbodycontours,zftblobs& ptFishblobs);

int getFishBlobCentreAndOrientation(cv::Mat imgFishAnterior,cv::Point2f ptCentre,int Angle,cv::Point2f& ptRevised,int& RevisedAngle);
int findMatchingContour(std::vector<std::vector<cv::Point> >& contours,
                              std::vector<cv::Vec4i>& hierarchy,
                              cv::Point pt,
                              int level);

#endif // FGMASKPROCESSING_H
