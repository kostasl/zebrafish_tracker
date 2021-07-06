/// \brief This class will contain the learning and heurestic required to detect fish and manage the update of their respective model instances
///
#include "fishdetector.h"
#include "template_detect.h"

static cv::Mat loadImage(const std::string& name)
{
    cv::Mat image = cv::imread(name, cv::IMREAD_GRAYSCALE);
    if (image.empty())
    {
        std::cerr << "Can't load image - " << name << std::endl;
        exit(-1);
    }
    return image;
}


fishdetector::fishdetector()
{
    mW_L1 = loadImage(std::string("/home/kostasl/workspace/zebrafishtrack/Rplots/KC_SparseNet.pgm"));
    mW_L2 = loadImage(std::string("/home/kostasl/workspace/zebrafishtrack/Rplots/outputLayer_trained.pgm") );

    cv::threshold(mW_L1,mW_L1,0.1,1,cv::THRESH_BINARY);
    cv::threshold(mW_L2,mW_L2,0.1,1,cv::THRESH_BINARY);

    mW_L1.convertTo(mW_L1, CV_32FC1);
    mW_L2.convertTo(mW_L2, CV_32FC1);
}


float fishdetector::netDetect(cv::Mat imgRegion)
{

    cv::Mat imgRegion_bin;
    cv::adaptiveThreshold(imgRegion,imgRegion_bin,1,cv::ADAPTIVE_THRESH_MEAN_C,cv::THRESH_BINARY,5,0);
    cv::imshow("input Bin",imgRegion_bin*255);

    // Input Is converted to Row Vector So we can do Matrix Multiplation
    assert(imgRegion_bin.cols*imgRegion_bin.rows == mW_L1.rows);
    cv::Mat vIn = imgRegion_bin.reshape(1,mW_L1.rows).t();  //Row Vector
    vIn.convertTo(vIn, CV_32FC1);

    // operation multiplies matrix A of size [a x b] with matrix B of size [b x c]
    //to the Layer 1 output produce matrix C of size [a x c]
    cv::Mat mL1_out = vIn*mW_L1;
    double minL1,maxL1;
    cv::minMaxLoc(mL1_out,&minL1,&maxL1);

    // Threshold for Activation Function
    cv::threshold(mL1_out,mL1_out,fL1_activity_thres,1,cv::THRESH_BINARY);


    cv::Mat mL2_out =  mL1_out*mW_L2.t();

    //cv::imshow("L1 Out", mL1_out);
    //Output fraction of Active Input that is filtered by Synaptic Weights, (Fraction of Active Pass-through KC neurons)
    float fOut = mL2_out.at<float>(0,0)/mW_L1.cols;
    qDebug() << "***R: " << fOut << " KCmin: "<< minL1 << " KCmax: " << maxL1;
    //cv::imshow("output",mL2_out);

    return(fOut);
}

