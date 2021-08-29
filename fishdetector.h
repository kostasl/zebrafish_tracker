#ifndef FISHDETECTOR_H
#define FISHDETECTOR_H

//#include <config.h>
//#include <fishmodel.h>
//#include <larvatrack.h>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
#include <opencv2/highgui/highgui.hpp>



typedef cv::KeyPoint zftblob;
typedef std::vector<zftblob> zftblobs;

/// /brief An Random Kernel - Mushroom body inspiried - recognition neural net implementation
/// Matrices are learned offline in R script that uses Image Examples to set NN output layer weights
/// see tracker_img_recognitionNN.R
class fishdetector
{
public:
    fishdetector();
    float netNeuralTF(float a);
    static cv::Mat getNormedBoundedImg(cv::Mat& frame, cv::RotatedRect fishRotAnteriorBox); //Normed Bounded region of Rotated Rect
    static cv::Mat getNormedTemplateImg(cv::Mat& frame, cv::RotatedRect& fishRotAnteriorBox); // Normed Rot Rect Image
    float netDetect(cv::Mat imgRegion_bin,float &fFishClass,float & fNonFishClass);
    float scoreBlobRegion(cv::Mat frame,zftblob& fishblob,cv::Mat& outframeAnterior_Norm,cv::Mat& outmaskRegionScore,std::string regTag);
    float fL1_activity_thres = 10; //# Number of INput that need to be active for KC to fire/Activate

private:

    cv::Mat mW_L1,mB_L1; //Weights and Biases Layer 1
    cv::Mat mW_L2,mB_L2; //Weights and Biases Layer 2
    cv::Mat mL1_out; ///Matrix Holding Result of Input*L1
    cv::Mat mL2_out;
};

#endif // FISHDETECTOR_H
