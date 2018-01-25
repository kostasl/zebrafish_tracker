#ifndef FGMASKPROCESSING_H
#define FGMASKPROCESSING_H

#include <opencv2/opencv.hpp>
#include <GUI/mainwindow.h>
#include <QString>
#include <QApplication>

unsigned int getBGModelFromVideo(cv::Mat& fgMask,MainWindow& window_main,QString videoFilename,QString outFileCSV,unsigned int MOGhistoryLength);
bool updateBGFrame(cv::Mat& frame, cv::Mat& fgMask, unsigned int nFrame,uint MOGhistory);

#endif // FGMASKPROCESSING_H
