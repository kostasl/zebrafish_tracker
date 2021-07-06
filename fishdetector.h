#ifndef FISHDETECTOR_H
#define FISHDETECTOR_H

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
#include <opencv2/highgui/highgui.hpp>

/// /brief An Random Kernel - Mushroom body inspiried - recognition neural net implementation
/// Matrices are learned offline in R script that uses Image Examples to set NN output layer weights
/// see tracker_img_recognitionNN.R
class fishdetector
{
public:
    fishdetector();
    float netDetect(cv::Mat imgRegion);
    float fL1_activity_thres = 10; //# Number of INput that need to be active for KC to fire/Activate

private:
    cv::Mat mW_L1;
    cv::Mat mW_L2;
};

#endif // FISHDETECTOR_H
