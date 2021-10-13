#ifndef FISHDETECTOR_H
#define FISHDETECTOR_H

//#include <config.h>
//#include <fishmodel.h>
//#include <larvatrack.h>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
#include <opencv2/highgui/highgui.hpp>

#include <tensorDNN/tf_image.hpp>

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
    static cv::Mat getNormedBoundedImg(const cv::Mat& frame, cv::RotatedRect fishRotAnteriorBox); //Normed Bounded region of Rotated Rect
    static cv::Mat getNormedTemplateImg(const cv::Mat& frame, cv::RotatedRect& fishRotAnteriorBox); // Normed Rot Rect Image
    float netDetect(cv::Mat imgRegion_bin,float &fFishClass,float & fNonFishClass);
    float netDNNDetect(cv::Mat imgRegion_bin,float &fFishClass,float & fNonFishClass);

    float scoreBlobRegion(cv::Mat frame,zftblob& fishblob,cv::Mat& outframeAnterior_Norm,cv::Mat& outmaskRegionScore,std::string regTag);
    float fL1_activity_thres = 10; //# Number of INput that need to be active for KC to fire/Activate
    static void test(); //Test with given Images
    static void testTFModelPrediction(const std::vector<cv::Mat>& image); //Test  Tensorflow DNN model on provided Image
private:

    /// Vector of Matrices Holding Layer state / Result of FWD prop
    //Weights and Biases Layer
    std::vector<cv::Mat> vmW_L;
    std::vector<cv::Mat> vmB_L;
    std::vector<cv::Mat> vmL_out;

    // TensorFlow Assistant Class -
    // Only 20% of the available GPU memory will be allocated
    float m_gpu_memory_fraction = 0.2f;
    // the model will try to infer the input and output layer names automatically
    // (only use if it's a simple "one-input -> one-output" model
    bool m_inferInputOutput = false;
    const std::string mSavedModelPath = "/home/kostasl/workspace/zebrafishtrack/tensorDNN/savedmodels/fishNet_prob/";
    tf_image::TF_Model m_TFmodel;



    //cv::Mat mL2_out;
    //cv::Mat mL3_out;
    //cv::Mat mL4_out;
    //cv::Mat mL5_out;
};

#endif // FISHDETECTOR_H
