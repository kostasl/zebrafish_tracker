#ifndef TEMPLATE_DETECT_H
#define TEMPLATE_DETECT_H


#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
#include <opencv2/highgui/highgui.hpp>
#include <QString>



///
/// \brief loadTemplatesFromDirectory
/// \param strDir
/// \return
///
int loadTemplatesFromDirectory(QString strDir);

//makeTemplateCache(angleStep)
void makeTemplateVar(cv::Mat& templateIn,cv::Mat& imgTemplateOut, int iAngleStepDeg);
//int templatefindFishInImage(cv::Mat& imgGreyIn,cv::Mat& imgtempl,cv::Size templSz, double& matchScore, cv::Point& locations_tl,int& startRow,int& startCol);
int templatefindFishInImage(cv::UMat& imgGreyIn,cv::UMat& imgtempl,cv::Size templSz, double& matchScore, cv::Point& locations_tl,int& startRow,int& startCol,bool findFirstMatch);
double doTemplateMatchAroundPoint(const cv::Mat& maskedImg_gray,cv::Point pt,int& detectedAngle,cv::Point& detectedPoint ,cv::Mat& frameOut );

///
/// \brief addTemplateToCache Expands the Cache with all the fish body image templates by one row
///        Assumes all template images have the same size
/// \param imgTempl
/// \param idxTempl
/// \return idxTempl
///
int addTemplateToCache(cv::Mat& imgTempl,cv::UMat& FishTemplateCache,int idxTempl);


int deleteTemplateRow(cv::Mat& imgTempl,cv::UMat& FishTemplateCache,int idxTempl);
//expandTemplateCache(newTemplateImage)
//templatefindFishInImage()
#endif // TEMPLATE_DETECT_H
