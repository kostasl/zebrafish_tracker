#ifndef TEMPLATE_DETECT_H
#define TEMPLATE_DETECT_H


#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
#include <opencv2/highgui/highgui.hpp>




//makeTemplateCache(angleStep)
void makeTemplateVar(cv::Mat& templateIn,cv::Mat& imgTemplateOut, int iAngleStepDeg);
int templatefindFishInImage(cv::Mat& imgGreyIn, cv::Mat& imgtempl, cv::Size templSz, double& matchScore, cv::Point& locations_tl, int& startRow);

//expandTemplateCache(newTemplateImage)
//templatefindFishInImage()
#endif // TEMPLATE_DETECT_H
