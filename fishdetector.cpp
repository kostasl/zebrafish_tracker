/// \brief Class utilizing trained neural network that detects fish anterior at size of a template img and
/// corrects the position the blob centre to so as to assist correct eye detection.
/// The NN classifier is trained in tracker_img_recognitionNN.R script using a collection of fish and non-fish images and saved as Matrices
/// exported as YAML in fishNet.yml. These are loaded as OpenCV matrices and is used to classify candidate fish blobs.
///
///
///

//#include <config.h>
//#include "larvatrack.h"

//#include "fishmodel.h"


#include "fishdetector.h"
#include "template_detect.h"

#include <tensorflow/c/c_api.h> // TensorFlow C API header.
#include <tensorDNN/tf_image.hpp>

extern trackerState gTrackerState;
extern cv::Mat gFishTemplateCache; //A mosaic image contaning copies of template across different angles

//TF Aux funct
void NoOpDeallocator(void* data, size_t a, void* b) {}


static cv::Point2f rotate2d(const cv::Point2f& inPoint, const double& angRad)
{
    cv::Point2f outPoint;
    //CW rotation
    outPoint.x = std::cos(angRad)*inPoint.x - std::sin(angRad)*inPoint.y;
    outPoint.y = std::sin(angRad)*inPoint.x + std::cos(angRad)*inPoint.y;
    return outPoint;
}


static cv::Point2f rotateAboutPoint(const cv::Point2f& inPoint, const cv::Point2f& center, const double& angRad)
{
    return rotate2d(inPoint - center, angRad) + center;
}


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
    String sDir = std::string("/home/kostasl/workspace/zebrafishtrack/Rplots/fishNet.yml");

    FileStorage fsNet;
    fsNet.open( sDir, FileStorage::READ,String("UTF-8"));
    int iLayerCount = 0;

    // Load Network Into Array of Matrices -
    fsNet["NLayer"] >> iLayerCount;
    vmW_L.resize(iLayerCount);
    vmB_L.resize(iLayerCount);
    vmL_out.resize(iLayerCount);

    ///Matrices Are read Serialized per Column - not per row (as is the default R ser. of matrices)
    for (int l=0;l<iLayerCount;l++)
    {
        //cv::Mat mW_L;
        //cv::Mat mB_L;
        string sLayerW_ID = string("LW") + QString::number(l+1).toStdString();
        string sLayerB_ID = string("LB") + QString::number(l+1).toStdString();

        fsNet[sLayerW_ID] >> vmW_L[l];
        fsNet[sLayerB_ID] >> vmB_L[l];

        //vmW_L[l].convertTo(vmW_L[l], CV_32FC1);
        //vmB_L[l].convertTo(vmW_L[l], CV_32FC1);
    }

    //fsNet["LW2"] >> mW_L2;    fsNet["LW3"] >> mW_L3;    fsNet["LW4"] >> mW_L4;    fsNet["LW5"] >> mW_L5;    fsNet["LB1"] >> mB_L1;    fsNet["LB2"] >> mB_L2;    fsNet["LB3"] >> mB_L3;
    //fsNet["LB4"] >> mB_L4;
    //fsNet["LB5"] >> mB_L5;

    //mW_L1.convertTo(mW_L1, CV_32FC1);    mW_L2.convertTo(mW_L2, CV_32FC1);    mW_L3.convertTo(mW_L3, CV_32FC1);    mW_L4.convertTo(mW_L4, CV_32FC1);
    //mB_L1.convertTo(mB_L1, CV_32FC1);    mB_L2.convertTo(mB_L2, CV_32FC1);    mB_L3.convertTo(mB_L3, CV_32FC1);    mB_L4.convertTo(mB_L4, CV_32FC1);
}

/// \brief Utility function takes img contained in rotated rect and returns the contained image region
/// rotated up-right (orthonormal ) - Used to obtain templates to train classifier
cv::Mat fishdetector::getNormedBoundedImg(const cv::Mat& frame, cv::RotatedRect fishRotAnteriorBox)
{
    cv::Mat imgFishAnterior, imgFishAnterior_Norm;
    /// Size Of Norm Head Image
    cv::Rect fishBoundingRect = fishRotAnteriorBox.boundingRect();
    // fishBoundingRect.width  +=2; fishBoundingRect.height +=2;
    cv::Size szFishAnteriorNorm = cv::Size(max (fishBoundingRect.size().width,fishBoundingRect.size().height),
                                           max (fishBoundingRect.size().width,fishBoundingRect.size().height));//

    /// Extract Region and rotate to orient template region vertically

    //Threshold The Match Check Bounds Within Image
    cv::Rect imgBounds(0,0,frame.cols,frame.rows);


    if (!( //Check IMage Bounds contain the whole bounding box//
        imgBounds.contains(fishBoundingRect.br()) &&
            imgBounds.contains(fishBoundingRect.tl())))
        return (imgFishAnterior_Norm); //This region Is out Of Bounds /

    // Use the FG Image to extract Head Frame
    frame(fishBoundingRect).copyTo(imgFishAnterior);

    /// Make Rotation MAtrix About Centre Of Cropped Image
    cv::Point2f ptRotCenter = fishRotAnteriorBox.center - fishRotAnteriorBox.boundingRect2f().tl();

    cv::Mat Mrot = cv::getRotationMatrix2D( ptRotCenter, fishRotAnteriorBox.angle,1.0); //Rotate Upwards

    ///Make Rotation Transformation
    //Need to fix size of Upright/Normed Image
    cv::warpAffine(imgFishAnterior,imgFishAnterior_Norm,Mrot,szFishAnteriorNorm);

    //Make Sure Normed Template Fits in Bounded Region
    //assert(imgFishAnterior_Norm.cols >= fishRotAnteriorBox.size.width);
    //assert(imgFishAnterior_Norm.rows >= fishRotAnteriorBox.size.height);

    return imgFishAnterior_Norm;
}

/// Extracts a template sized image region contained in rotatedRect and returns the image Vert orientated - Normalized Template IMg
cv::Mat fishdetector::getNormedTemplateImg(const cv::Mat& frame, cv::RotatedRect& fishRotAnteriorBox)
{
    cv::Size szFishAnteriorNorm = fishRotAnteriorBox.boundingRect().size();
    cv::Mat imgFishAnterior_Norm;


    //Define Regions and Sizes for extracting Orthonormal Fish
    //Top Left Corner of templateSized Rect relative to Rectangle Centered in Normed Img
    cv::Size szTemplateImg = gTrackerState.gszTemplateImg;

    cv::Point ptTopLeftTemplate(min(szFishAnteriorNorm.width, max(0,szFishAnteriorNorm.width/2-szTemplateImg.width/2)),
                                min(szFishAnteriorNorm.height, max(0,szFishAnteriorNorm.height/2-szTemplateImg.height/2)) );

    cv::Rect rectFishTemplateBound = cv::Rect(ptTopLeftTemplate,szTemplateImg);

    cv::Mat imgBoundedNorm = getNormedBoundedImg(frame, fishRotAnteriorBox);

    // RETURN EMPTY If FishAnterior Image is too small (Near boundary case)
    if ((rectFishTemplateBound.width + rectFishTemplateBound.x) > imgBoundedNorm.cols)
        return imgFishAnterior_Norm;
    if ((rectFishTemplateBound.height + rectFishTemplateBound.y) > imgBoundedNorm.rows)
        return imgFishAnterior_Norm;


    ///Cut Down To Template Size - Take The sized Template Region at the centre of the rotated image
    imgFishAnterior_Norm = imgBoundedNorm(rectFishTemplateBound);

    return(imgFishAnterior_Norm);
}


///brief
///
///

cv::Mat sparseBinarize(cv::Mat& imgRegion,float targetdensity)
{
    cv::Mat imgRegion_bin,imgRegion_bin_dense;

    /// Sparse Binarize Through Adaptive Threshold to enhance fish-like pattern // Substract - const val from mean
    // Ideally We want to maintain input sparseness
    double activePixRatio = 1.0;
    int meanAdapt = 0;
    int maxIter = 20;

      /// Regulate Input Sparseness
        while (activePixRatio >  targetdensity &&
               maxIter > 0)
        {
          cv::adaptiveThreshold(imgRegion,imgRegion_bin,1,cv::ADAPTIVE_THRESH_MEAN_C,cv::THRESH_BINARY,imgRegion.size().width,meanAdapt);
          activePixRatio = (cv::sum(imgRegion_bin)[0])/(imgRegion_bin.cols*imgRegion_bin.rows);
          meanAdapt += (int)255*(targetdensity-activePixRatio)-2;
          if (maxIter == 20) //Save First Binarized most Dense
              imgRegion_bin.copyTo(imgRegion_bin_dense);
          maxIter--;
        }

    return(imgRegion_bin);
}

/// \brief Two step classificiation of region : First, it uses Neural
/// Net to scan region around provided blog and provide a detection score as a mask (returns max score value)
/// If blob passes the FishNet classification threshold , the blobs centre position is changed to the point of max classification score.
/// \todo could do image Pyramids to scan Across Scales
/// @param regTag an Id for debugging purposes
/// @outframeAnterior_Norm returns image of isolated head centered at best detection point according to NN

float fishdetector::scoreBlobRegion(cv::Mat frame,zftblob& fishblob,cv::Mat& outframeAnterior_Norm,
                                    cv::Mat& outmaskRegionScore,string regTag="0")
{
  cv::Mat imgFishAnterior,imgFishAnterior_Norm,imgFishAnterior_Norm_bin,imgFishAnterior_Norm_tmplcrop;
  cv::Mat imgFishAnterior_Norm_bin_dense; //Used to Find Contours

  if ((fishblob.pt.x + gTrackerState.gFishBoundBoxSize) > frame.cols ||
      (fishblob.pt.x - gTrackerState.gFishBoundBoxSize) < 0)
      return(0.0f);
  if ((fishblob.pt.y + gTrackerState.gFishBoundBoxSize) > frame.rows ||
          (fishblob.pt.y - gTrackerState.gFishBoundBoxSize) < 0 )
      return(0.0f);

  // Take bounding of blob,
  // get Rotated Box Centre Coords relative to the cut-out of the anterior Body - This we use to rotate the image
  cv::RotatedRect fishRotAnteriorBox(fishblob.pt,
                                      cv::Size(gTrackerState.gFishBoundBoxSize,gTrackerState.gFishBoundBoxSize),
                                                fishblob.angle);
   /// Size Of Norm Head Image
   cv::Size szFishAnteriorNorm = fishRotAnteriorBox.boundingRect().size();// (min(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height)+4,                              max(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height)+4);
   // Define SCore Canvas - Where we draw results from Recognition Scanning
   cv::Mat maskRegionScore_Norm((int)szFishAnteriorNorm.width, (int)szFishAnteriorNorm.height, CV_32FC1, cv::Scalar(0, 0, 0));


   //To Check Bounds Within Image
   cv::Rect imgBounds(0,0,frame.cols,frame.rows);

   // Check if region size is large enough to scan for fish
   if (!( //Looks Like a fish is found, now Check Bounds
       imgBounds.contains(fishRotAnteriorBox.boundingRect().br()) &&
           imgBounds.contains(fishRotAnteriorBox.boundingRect().tl())))
            return (-1); //This Fish Is out Of Bounds /


   imgFishAnterior_Norm =  getNormedBoundedImg(frame,fishRotAnteriorBox);
   //cv::imshow("Test Tmpl Norm",imgTemplate_Norm);

  /// Extract Region and rotate to orient larva body vertically
  // Use the FG Image to extract Head Frame
  //frame(fishRotAnteriorBox.boundingRect()).copyTo(imgFishAnterior);

  cv::Point2f ptRotCenter = fishRotAnteriorBox.center - fishRotAnteriorBox.boundingRect2f().tl();

  //Binarize Input To set Specific Sparseness/Density
  //imgFishAnterior_Norm_bin = sparseBinarize(imgFishAnterior_Norm,gTrackerState.fishnet_inputSparseness);
  cv::normalize(imgFishAnterior_Norm,imgFishAnterior_Norm_bin,1.0,0,NORM_MINMAX,CV_32FC1);
  //cv::imshow(std::string("BIN_N") + regTag,imgFishAnterior_Norm_bin);


  /// SliDing Window Scanning
  int iSlidePx_H_step = 2;
  int iSlidePx_H_begin = ptRotCenter.x- gTrackerState.gszTemplateImg.width/2 - 4;//max(0, imgFishAnterior_Norm.cols/2- sztemplate.width);
  int iSlidePx_H_lim = iSlidePx_H_begin+4;  //imgFishAnterior_Norm.cols/2; //min(imgFishAnterior_Norm.cols-sztemplate.width, max(0,imgFishAnterior_Norm.cols/2+ sztemplate.width) ) ;

   // V step - scanning for fishhead like image in steps
  int iSlidePx_V_step = 2;
  int iSlidePx_V_begin = std::max(0,(int)(ptRotCenter.y - gTrackerState.gszTemplateImg.height/2)-8); //(int)(ptRotCenter.y - sztemplate.height) sztemplate.height/2
  int iSlidePx_V_lim = iSlidePx_V_begin + 8;//min(imgFishAnterior_Norm.rows - gTrackerState.gszTemplateImg.height, iSlidePx_V_begin + 10); //(int)(sztemplate.height/2)


  float scoreFish,scoreNonFish,dscore; //Recognition Score tested in both Vertical Directions
  // Do netDetect using a Sliding window
  cv::Mat imgFishAnterior_Norm_tmplcrop_vflip;
  for (int i=iSlidePx_H_begin;i <= iSlidePx_H_lim;i+=iSlidePx_H_step)
  {
      for (int j=iSlidePx_V_begin;j <= iSlidePx_V_lim;j+=iSlidePx_V_step)
      {
          //Define Regions Top-Left Corner Position within Extracted Blob Region and
          cv::Point ptTopLeftTemplate(i,j);
          cv::Rect rectFishTemplateBound = cv::Rect(ptTopLeftTemplate, gTrackerState.gszTemplateImg );
          //CROP Extract a Template sized subregion of Orthonormal Fish
          imgFishAnterior_Norm_bin(rectFishTemplateBound).copyTo(imgFishAnterior_Norm_tmplcrop);

          dscore = this->netDetect(imgFishAnterior_Norm_tmplcrop,scoreFish,scoreNonFish);

          ///Do Not Test Orientation - Blob Should Have the correct Angle
          //Check Both Vertical Orientations
          //cv::flip(imgFishAnterior_Norm_tmplcrop, imgFishAnterior_Norm_tmplcrop_vflip, 0);
          /// Score Number of Bin Pix In cropped region so as to push for regions that fully contain the eyes
          //float activePixRatio = (1+cv::sum(imgFishAnterior_Norm_tmplcrop)[0])/(imgFishAnterior_Norm_tmplcrop.cols*imgFishAnterior_Norm_tmplcrop.rows);
          ///  Store recognition score in Mask at(row,col) -//
          maskRegionScore_Norm.at<float>(j +gTrackerState.gszTemplateImg.height/2, // ,
                                         i +gTrackerState.gszTemplateImg.width/2) = dscore;// //max(0.0f,dscore);//(scoreFish + scoreNonFish + 1e-3))/activePixRatio; //+ activePixRatio; //
          //qDebug() << "(" << i+sztemplate.width/2 << "," <<j+sztemplate.height/2<<") = " << round(sc1*100)/100.0;
            //qDebug() << maskRegionScore_Norm.at<float>(j+gTrackerState.gszTemplateImg.height/2,
            //                                           i+gTrackerState.gszTemplateImg.width/2);
        }//For Each Vertical
    }//For Each Horizontal

    /// Find and Best scorring point ScoreMask - Along with Normed Fish Image centered at best-point
    //Need to Rotate Score Image Back to Original Orientiation to return Coordinates oF Best Match
    /// Find Max Match Position In Non-Norm pic (original orientation)
    double minL1,maxL1;
    cv::Point ptmin,ptmax;
    cv::GaussianBlur(maskRegionScore_Norm,maskRegionScore_Norm,cv::Size(3,3),3,3);
    cv::minMaxLoc(maskRegionScore_Norm,&minL1,&maxL1,&ptmin,&ptmax);
    // Rotate Max Point Back to Original Orientation
    //cv::Mat MrotInv = cv::getRotationMatrix2D( ptRotCenter, -fishblob.angle,1.0); //Rotate Upwarte Upwards
    //cv::warpAffine(maskRegionScore,outmaskRegionScore,MrotInv,szFishAnteriorNorm);
    cv::Point ptmax_orig = rotateAboutPoint(ptmax,cv::Point2f(maskRegionScore_Norm.cols/2,maskRegionScore_Norm.rows/2),
                                            (fishblob.angle)*(CV_PI/180.0) ); //-fishblob.angle   angle
    // End Of FishNet Detection //


    //Update Blob Location And add Classifier Score
    fishblob.response = maxL1; //Save Recognition Score
    if (maxL1 > 0)
        fishblob.pt = ptmax_orig+fishRotAnteriorBox.boundingRect().tl(); //Shift Blob Position To Max  To Max Recognition Point

    // DEBUG IMG //
    //cv::circle(imgFishAnterior,ptmax_orig,4,CV_RGB(250,200,210),2);
    //cv::imshow(string("Fish Region Body ") + regTag,imgFishAnterior);
    // DEBUG IMG //
    //cv::circle(imgFishAnterior_Norm,ptmax,3,CV_RGB(200,200,210),2);
    //cv::circle(maskRegionScore_Norm,ptmax,3,CV_RGB(0,0,0),2);


    cv::imshow(string("Fish Region Body Norm ") + regTag,imgFishAnterior_Norm_bin);
    // DEBUG IMG //
    //cv::normalize(maskRegionScore_Norm, maskRegionScore_Norm, 0, 1, cv::NORM_MINMAX);
    //cv::imshow(string("Score Mask Body Norm") + regTag,maskRegionScore_Norm);                                   gionScore_Norm);

    /// Find Max Score Coords In Normed FishAnterior / Around Best Match Region (Using Normed Regiest Match Region (Using Normed Region)

    cv::Point ptTopLeftTemplate(max(0,ptmax.x-gTrackerState.gszTemplateImg.width/2),      //   TemplateImg.width/2),
                                max(0,ptmax.y-gTrackerState.gszTemplateImg.height/2) );
    // Stick To Boundary for Template Size Window
    ptTopLeftTemplate.x =((ptTopLeftTemplate.x + gTrackerState.gszTemplateImg.width) >= maskRegionScore_Norm.cols)?maskRegionScore_Norm.cols-gTrackerState.gszTemplateImg.width:ptTopLeftTemplate.x;
    ptTopLeftTemplate.y =((ptTopLeftTemplate.y + gTrackerState.gszTemplateImg.height) >= maskRegionScore_Norm.rows)?maskRegionScore_Norm.rows-gTrackerState.gszTemplateImg.height: ptTopLeftTemplate.y;
    // OR Simply Take Midline Of Image as best Guess of Fish Location
    //cv::Point ptTopLeftTemplate(szFishAnteriorNorm.width/2-gTrackerState.gLastfishimg_template.size().width/2,
    //                         szFishAnteriorNorm.height/2-gTrackerState.gLastfishimg_template.size().height/2);

//    // Bound search region so cropping Of Norm does not exceed image bounds
//    cv::Rect rect_bound = cv::Rect(cv::Point(gTrackerState.gLastfishimg_template.size().width/2,
//                                             gTrackerState.gLastfishimg_template.size().height/2),
//                                             cv::Size(maskRegionScore.cols-gTrackerState.gLastfishimg_template.size().width/2,
//                                                       maskRegionScore.rows-gTrackerState.gLastfishimg_template.size().height/2));
//    cv::minMaxLoc(maskRegionScore(rect_bound),&minL1,&maxL1,&ptmin,&ptmax);
//    cv::Point ptTopLeftTemplate(min(imgFishAnterior_Norm.cols, max(0, ptmax.x - gTrackerState.gLastfishimg_template.size().width/2)),
//                                min(imgFishAnterior_Norm.rows, max(0, ptmax.y - gTrackerState.gLastfishimg_template.size().height/2))
//                               );


    cv::Rect rectFishTemplateBound = cv::Rect(ptTopLeftTemplate, gTrackerState.gszTemplateImg );
    assert(rectFishTemplateBound.width + rectFishTemplateBound.x <= imgFishAnterior_Norm.cols );
    assert(rectFishTemplateBound.height + rectFishTemplateBound.y <= imgFishAnterior_Norm.rows );

    /// CROP Extract a Template sized subregion of Orthonormal Fish ///
    imgFishAnterior_Norm_tmplcrop       = imgFishAnterior_Norm(rectFishTemplateBound);
    imgFishAnterior_Norm_tmplcrop.copyTo(outframeAnterior_Norm);
    maskRegionScore_Norm(rectFishTemplateBound).copyTo(outmaskRegionScore);

    /// Set Mark Point For Eye Detection ///
    //if (gTrackerState.bAdaptEyeMaskVOffset)
    //    gTrackerState.giHeadIsolationMaskVOffset = min(maxpt.y,imgFishAnterior_Norm.rows);//ptmax.y



  return (maxL1);

}

float fishdetector::netNeuralTF(float a)
{
    return( (1.0f/( 1.0f+ std::exp(-a) )) );
}
/// \brief Applies pre-trained MB like NN on Binarized Input image
/// Networks supports two L2 neurons - These recognition nets suffer from decreases in input sparseness :
/// More active inputs increase the output scores - reducing the networks selectivity
float fishdetector::netDetect(cv::Mat imgRegion_bin,float &fFishClass,float & fNonFishClass)
{
    //fL1_activity_thres = gTrackerState.fishnet_L1_threshold;

    // Input Is converted to Row Vector So we can do Matrix Multiplation
    //assert(imgRegion_bin.cols*imgRegion_bin.rows == mW_L1.cols);
    if (imgRegion_bin.cols*imgRegion_bin.rows != vmW_L[0].cols)
    {
        fFishClass = -1;
        fNonFishClass = -1;
        return 0.0;
    }

    //qDebug() << "input img Type:" << type2str(imgRegion_bin.type());
    //cv::imshow("Input Img ", imgRegion_bin);

    // Make Col Vector
    cv::Mat vIn(vmW_L[0].cols,1,CV_32FC1); /// =imgRegion_bin.reshape(0,mW_L1.cols);  <- Reshape Does not work as intented

    //for (int i=0; i<vIn.rows;i++){
    //      qDebug() << i << ". " << vIn.at<uchar>(0,i);
    //}
     // Reshape Image Matrix into column vector  (Custom required)
     int i = 0;
     for (int c=0; c<imgRegion_bin.cols;c++){
         for (int r=0; r<imgRegion_bin.rows;r++){
        vIn.at<float>(i,0) = (float)imgRegion_bin.at<uchar>(r,c)/255.0f;//netNeuralTF(vIn.at<float>(0,i));
        //qDebug() << i << ". " << vIn.at<float>(i,0);
        i++;
        }
    }
    //cv::imshow("Vin Out 8bit ", vIn.reshape(1,imgRegion_bin.rows));

    //vIn.convertTo(vIn, CV_32FC1);

    for (int l=0;l<vmW_L.size();l++)
    {
        if (l==0)
            vmL_out[l] = vmW_L[l]*vIn +vmB_L[l];
        else
            vmL_out[l] = vmW_L[l]*vmL_out[l-1] +vmB_L[l];

        //Apply Neural Transfer Function
        for (int i=0; i<vmL_out[l].rows;i++)
        {
            vmL_out[l].at<float>(i,0) = netNeuralTF(vmL_out[l].at<float>(i,0));
        }
    }
    //cv::imshow("Vin Out 32Fbit ", vIn.reshape(1,imgRegion_bin.rows));
    //qDebug() << "W_l1 Type:" << type2str(mW_L1.type());
    //cv::imshow("Vin Out TF", vIn.reshape(1,imgRegion_bin.rows));

//    /// \TODO Matrices are not read correctly beyond 1st column mW_L1
//    // operation multiplies matrix A of size [a x b] with matrix B of size [b x c]
//    //to the Layer 1 output produce matrix C of size [a x c]
//    mL1_out = mW_L1*vIn + mB_L1;
//    //Apply Neural Transfer Function
//    for (int i=0; i<mL1_out.cols;i++)
//    {
//        mL1_out.at<float>(0,i) = netNeuralTF(mL1_out.at<float>(0,i));
//    }

//    //Calc Layer 2 (Hidden Layer 2) Activation
//    mL2_out =  mW_L2*mL1_out + mB_L2;
//    //Apply Neural Transfer Function
//    for (int i=0; i<mL2_out.rows;i++)
//        mL2_out.at<float>(i,0) = netNeuralTF(mL2_out.at<float>(i,0));

//    //Calc Layer 3 (Output) Activation
//    mL3_out =  mW_L3*mL2_out + mB_L3;
//    //Apply Neural Transfer Function
//    for (int i=0; i<mL3_out.rows;i++)
//        mL3_out.at<float>(i,0) = netNeuralTF(mL3_out.at<float>(i,0));

//    //Calc Layer 4 (Output) Activation
//    mL4_out =  mW_L4*mL3_out + mB_L4;
//    //Apply Neural Transfer Function
//    for (int i=0; i<mL4_out.rows;i++)
//        mL4_out.at<float>(i,0) = netNeuralTF(mL4_out.at<float>(i,0));

//    //Calc Layer 5 (Output) Activation
//    mL5_out =  mW_L5*mL4_out + mB_L5;
//    //Apply Neural Transfer Function
//    for (int i=0; i<mL5_out.rows;i++)
//        mL5_out.at<float>(i,0) = netNeuralTF(mL5_out.at<float>(i,0));


//    //Output fraction of Active Input that is filtered by Synaptic Weights, (Fraction of Active Pass-through KC neurons)
//    fFishClass = mL5_out.at<float>(0,0);
//    // Check 2 row (neuron) output
//    fNonFishClass = mL5_out.at<float>(1,0);

    /// Net Output on Last layer
    fFishClass = vmL_out[vmL_out.size()-1].at<float>(0,0);
    fNonFishClass = vmL_out[vmL_out.size()-1].at<float>(1,0);

    //double minL1,maxL1;
    //cv::minMaxLoc(mL1_out,&minL1,&maxL1);
    //qDebug() << "***R: " << fOut << " KCmin: "<< minL1 << " KCmax: " << maxL1;
    //cv::imshow("output",mL2_out);

    return(fFishClass-fNonFishClass);
}

///\brief Unit test fishnet detector by loading and testing against the set of images used during training (R Script)
///  and comparing output to the networks output to the script in R
void fishdetector::test()
{
    qDebug() << "Testing fishNet Classifier Using stock template pics from disk.";
    QString strDirFish("/home/kostasl/workspace/zebrafishtrack/img/debug/fish/");
    QString strDirNonFish("/home/kostasl/workspace/zebrafishtrack/img/debug/nonfish/");

    std::vector<cv::Mat> vfish_mat = loadTemplatesFromDirectory(strDirFish);
    float fsumErrF =0.0f;
    float fsumErrNF =0.0f;
    float fCorrectF_class = 0.0;
    float fCorrectNF_class = 0.0;

    qDebug() << "~~~Test Fish templates~~~";
    float fishClassScore,nonfishScore,dscore;
    for (int i=0;i<vfish_mat.size();i++)
    {
        cv::Mat imgTempl = vfish_mat[i];

        dscore = gTrackerState.fishnet.netDetect(imgTempl,fishClassScore,nonfishScore);
        if (fishClassScore > -1)
            fsumErrF += pow((1-fishClassScore) + (0-nonfishScore),2);
        if (fishClassScore > nonfishScore) //Count Classification Errors
            fCorrectF_class += 1.0;
        else
           //  cv::imshow(string("Fail Fish")+to_string(i),imgTempl);
              cv::imwrite(gTrackerState.gstroutDirTemplates + string("/FailFish")+to_string(i)+".pgm", imgTempl);

         qDebug() << "Fish img gave F:" << fishClassScore << " NF:" << nonfishScore;
    }

    fsumErrF = fsumErrF/vfish_mat.size();
    fCorrectF_class = fCorrectF_class/vfish_mat.size();
    qDebug() << "Fish Class MSQ ERR:" << fsumErrF;

    qDebug() << "~~~Test NON-Fish templates~~~";
    std::vector<cv::Mat> vnonfish_mat = loadTemplatesFromDirectory(strDirNonFish);

    for (int i=0;i<vnonfish_mat.size();i++)
    {
        cv::Mat imgTempl = vnonfish_mat[i];

        dscore = gTrackerState.fishnet.netDetect(imgTempl,fishClassScore,nonfishScore);
        qDebug() << "Non-Fish img gave F:" << fishClassScore << " NF:" << nonfishScore;
        if (fishClassScore > -1)
            fsumErrNF += pow((0-fishClassScore) + (1-nonfishScore),2);
            if (nonfishScore >=fishClassScore)
            {
                fCorrectNF_class += 1.0;
            }else
               //cv::imshow(string("Fail NonFish")+to_string(i),imgTempl);
                cv::imwrite(gTrackerState.gstroutDirTemplates + string("/FailNONFish")+to_string(i)+".pgm", imgTempl);
    }
    fsumErrNF = fsumErrNF/vnonfish_mat.size();
    fCorrectNF_class = fCorrectNF_class/vnonfish_mat.size();

    qDebug() << "Non Fish Class MSQ ERR:" << fsumErrNF;

    qDebug() << "** Correctly classified %Fish:" << fCorrectF_class << " , %Non-Fish:" <<  fCorrectNF_class << " **";
    qDebug() << "~~~ Total MSQ Error:" << (fsumErrNF + fsumErrF)/2;

}


/// Example from https://github.com/kostasl/tensorflow_capi_sample
void fishdetector::testTFModelPrediction(const std::vector<cv::Mat>& vimages)
{
    // Only 20% of the available GPU memory will be allocated
    float gpu_memory_fraction = 0.2f;

    // the model will try to infer the input and output layer names automatically
    // (only use if it's a simple "one-input -> one-output" model
    bool inferInputOutput = false;

    // load a model from a .pb file
    tf_image::TF_Model model1;

    model1.loadModel("/home/kostasl/workspace/zebrafishtrack/tensorDNN/savedmodels/fishNet_prob/" , gpu_memory_fraction, inferInputOutput );//"graph_im2vec.pb"
    model1.setInputs( { "serving_default_sequential_1_input" } );
    model1.setOutputs( { "StatefulPartitionedCall" } );


    //cv::resize( image, image, { 28,38 } );

    // run prediction:
    for (int i=0; i< vimages.size();i++)
    {
        std::cout << i << ".  run prediction..." << std::endl;
        std::vector< std::vector< float > > results = model1.predict<std::vector<float>>( {vimages[i]} );
        //   ^              ^ second vector is a normal model output (i.e. for classification or regression)
        //   ^ the elements of the first vector correspond to the model's outputs (if the model has only one, the vector contains only 1 vector)

        // print results
        std::cout << "* print results n:" << results.size() << std::endl;
        for ( size_t i = 0; i < results.size(); i++ )
        {
          std::cout << "Output vector #" << i << ": ";
          for ( size_t j = 0; j < results[i].size(); j++ )
          {
            std::cout << std::fixed << std::setprecision(4) << results[i][j] << "\t";
          }
          if (results[i][0] > results[i][1])
              std::cout << "Image shows FISH"<< std::endl;
          else
              std::cout << "Image *NOT a FISH"<< std::endl;
          std::cout << std::endl;
        }
    }

}


