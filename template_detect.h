#ifndef TEMPLATE_DETECT_H
#define TEMPLATE_DETECT_H

#include <config.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
#include <opencv2/highgui/highgui.hpp>


/// CUDA //
#include <opencv2/opencv_modules.hpp> //THe Cuda Defines are in here
#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    #include "opencv2/cudaimgproc.hpp"
    #include "opencv2/cudaarithm.hpp"
    #include <opencv2/core/cuda.hpp>
#endif

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
int templatefindFishInImage(cv::Mat& imgGreyIn,cv::Mat& imgtempl,cv::Size templSz, double& matchScore, cv::Point& locations_tl,int& startRow,int& startCol,bool findFirstMatch);
double doTemplateMatchAroundPoint(const cv::Mat& maskedImg_gray,cv::Point pt,int& detectedAngle,cv::Point& detectedPoint ,cv::Mat& frameOut );

///
/// \brief addTemplateToCache Expands the Cache with all the fish body image templates by one row
///        Assumes all template images have the same size
/// \param imgTempl
/// \param idxTempl
/// \return idxTempl
///
int addTemplateToCache(cv::Mat& imgTempl,cv::Mat& FishTemplateCache,int idxTempl);


int deleteTemplateRow(cv::Mat& imgTempl,cv::Mat& FishTemplateCache,int idxTempl);

#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    double gpu_matchTemplate(cv::Mat templ_h,cv::cuda::GpuMat& image_h,cv::Point& ptBestMatch);
#endif

//expandTemplateCache(newTemplateImage)
//templatefindFishInImage()
#endif // TEMPLATE_DETECT_H
