/// \brief Class utilizing trained neural network that detects fish anterior at size of a template img and
/// corrects the position the blob centre to so as to assist correct eye detection.
///
/// \note
/// I tried two methods- Method A - A customly designed NN trained with BackProp : The NN classifier is trained in tracker_img_recognitionNN.R
/// script using a collection of fish and non-fish images and saved as Matrices
/// exported as YAML in fishNet.yml. These are loaded as OpenCV matrices and is used to classify candidate fish blobs.
/// The classifier perfomance of a 5 Layer version of this (or 7 layer) performace of this was poor when it came to the curved edges of the dish.
///
///   Method B: A DNN classifier using Tensorflow library : THis was trained using a python script found in tensorDNN subfolder.
///   The model is in the tensorDNN/savedmodels/fishNet and then it modified to add a softMax Output layer and saved again as tensorDNN/savedmodels/fishNet_prob
///   The tracker here uses an helper class TF_image (https://github.com/Xonxt/hello_tf_c_api) which I forked and had to modify
///   so it can load SavedModels without having to freeze them, thus utilizing the power of the newer TF V2.6 C API.
///
///
///

//#include <config.h>
#include "larvatrack.h"

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

/// Need to call before Any detection can occur - loads model files
bool fishdetector::initialize()
{

    String sDir = std::string(gTrackerState.strBackPropModelYMLFile);

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


    /// DNN tensorflow - There is py tool to get
    /// Load  Direction Agnostic model - used to detection position of fish in image region
    m_TFmodel_loc.loadModel( gTrackerState.strDNNTensorFlowModelFile, m_gpu_memory_fraction , m_inferInputOutput );//"graph_im2vec.pb"
    ///  You can obtain input/ouput names using saved_model_cli show --dir {mobilenet_save_path} --tag_set serve or by looking at the model name in python training script
    /// and then assume "serving_default_<model name>_input"
    m_TFmodel_loc.setInputs( { "serving_default_sequential_input" } );
    m_TFmodel_loc.setOutputs( { "StatefulPartitionedCall" } );

    ///\deprecated Load  Directional model - used to detection correct up-right image of fish in image region - used for direction detection
    //m_TFmodel_dir.loadModel(gTrackerState.strDNNTensorFlowVerticalModelFile ,m_gpu_memory_fraction , m_inferInputOutput );//"graph_im2vec.pb"
    //m_TFmodel_dir.setInputs( { "serving_default_sequential_1_input" } );
    //m_TFmodel_dir.setOutputs( { "StatefulPartitionedCall" } );
    bInitialized = true;
    return(true);
}


fishdetector::fishdetector()
{
    bInitialized = false;
}




/// \brief Utility function takes img contained in rotated rect and returns the contained image region
/// rotated up-right (orthonormal ) - Used to obtain templates to train classifier
/// Uses Contour To find and Correct Body Angle deviation IN Rotated Box
cv::Mat fishdetector::getNormedBoundedImg(const cv::Mat& frame, cv::RotatedRect fishRotAnteriorBox, bool correctOrientation = false)
{
    cv::Mat imgFishAnterior, imgFishAnterior_Norm,imgFishAnterior_bin;
    cv::Mat imgContourDBG;
    std::vector<std::vector<cv::Point> > fishAnteriorcontours;
    std::vector<cv::Vec4i> fishAnteriorhierarchy;
    cv::Rect fishBoundingRect = fishRotAnteriorBox.boundingRect();
    double rotationNormAngle;// = fishRotAnteriorBox.angle;
    cv::Point2f ptRotCenter;// = cv::Point2f(fishBoundingRect.tl()) - fishRotAnteriorBox.center;
    /// Size Of Norm Head Image

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

    if (correctOrientation)//If Correction Flag Is set
    {
        imgContourDBG = cv::Mat::zeros(imgFishAnterior.rows,imgFishAnterior.cols,imgFishAnterior.type());
        //Uses FG MASKED Image With Very High threshold To Find Swim Bladder Orientation
        cv::threshold(gTrackerState.mfgFishFrame(fishBoundingRect),imgFishAnterior_bin,max(30,gTrackerState.g_FGSegthresh*6),255,THRESH_BINARY);
        //cv::imshow("getNormedBoundedImg",imgFishAnterior_bin);
        cv::findContours( imgFishAnterior_bin, fishAnteriorcontours,fishAnteriorhierarchy, cv::RETR_CCOMP,
                          cv::CHAIN_APPROX_SIMPLE , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE
        int fishidx  = findMatchingContour(fishAnteriorcontours,fishAnteriorhierarchy,cv::Point(imgFishAnterior_bin.cols/2,imgFishAnterior_bin.rows/2),-1 );
        if (fishAnteriorcontours.size() > 0)
        {
            //drawContours(imgContourDBG,fishAnteriorcontours,fishidx,CV_RGB(255,255,255));
            //cv::imshow("getNormedBoundedImg_COntour",imgContourDBG);

            if (fishAnteriorcontours[fishidx].size() > 5)
            {
                cv::RotatedRect boundEllipse = cv::fitEllipse(fishAnteriorcontours[fishidx]);
                int DAngle = getAngleDiff(fishRotAnteriorBox.angle,boundEllipse.angle); //Correct Angle
                if (abs(DAngle) < 90){
                    fishRotAnteriorBox.angle += DAngle;
                    qDebug() << "Th.Fix:" << DAngle;
                }
            }
        }  //If Correction Flag Is set.. otherwise just copy to output
    }

    rotationNormAngle = fishRotAnteriorBox.angle;
    ptRotCenter = fishRotAnteriorBox.center - fishRotAnteriorBox.boundingRect2f().tl();

    /// Make Rotation MAtrix About Centre Of Cropped Image

    cv::Mat Mrot = cv::getRotationMatrix2D( ptRotCenter, rotationNormAngle,1.0); //Rotate Upwards
    ///Make Rotation Transformation
    //Need to fix size of Upright/Normed Image
    cv::warpAffine(imgFishAnterior,imgFishAnterior_Norm,Mrot,szFishAnteriorNorm);

    //Make Sure Normed Template Fits in Bounded Region
    //assert(imgFishAnterior_Norm.cols >= fishRotAnteriorBox.size.width);
    //assert(imgFishAnterior_Norm.rows >= fishRotAnteriorBox.size.height);

    return imgFishAnterior_Norm;
}

/// Extracts a template sized image region contained in rotatedRect and returns the image Vert orientated - Normalized Template IMg
cv::Mat fishdetector::getNormedTemplateImg(const cv::Mat& frame, cv::RotatedRect& fishRotAnteriorBox,bool correctOrientation = false)
{
    cv::Size szFishAnteriorNorm = fishRotAnteriorBox.boundingRect().size();
    cv::Mat imgFishAnterior_Norm;


    //Define Regions and Sizes for extracting Orthonormal Fish
    //Top Left Corner of templateSized Rect relative to Rectangle Centered in Normed Img
    cv::Size szTemplateImg = gTrackerState.gszTemplateImg;

    cv::Point ptTopLeftTemplate(min(szFishAnteriorNorm.width, max(0,szFishAnteriorNorm.width/2-szTemplateImg.width/2)),
                                min(szFishAnteriorNorm.height, max(0,szFishAnteriorNorm.height/2-szTemplateImg.height/2)) );

    cv::Rect rectFishTemplateBound = cv::Rect(ptTopLeftTemplate,szTemplateImg);

    cv::Mat imgBoundedNorm = getNormedBoundedImg(frame, fishRotAnteriorBox,correctOrientation);

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

/// \brief  classificiation of region : First, it uses Deep Neural Net
/// to scan region around provided blog and provide a detection score as a mask (returns max score value)
/// If blob passes the FishNet classification threshold , the blobs centre position is changed to the point of max classification score.
/// @param regTag an Id for debugging purposes
/// @param iSlidePx_H_step  Grid density Hz: Skip N px between scan points Horizontal
/// @param iSlidePx_V_step - Grid density V: Skip N px between scan points Vertically
/// @param iSlidepxLim - Number of px To scan in either direction (ie Scan Square Size)
/// @outframeAnterior_Norm returns image of isolated head centered at best detection point by the DNN and oriented according to blob angle
/// \note Only updates Blob position and response score if new classifier score is higher than existing
float fishdetector::scoreBlobRegion(cv::Mat frame,zftblob& fishblob,cv::RotatedRect fishRotAnteriorBox,cv::Mat& outframeAnterior_Norm,
                                    cv::Mat& outmaskRegionScore, int iSlidepxLim = 40,int iSlidePx_H_step = 10, int iSlidePx_V_step = 10, string regTag="0", bool bstopAtFirstMatch =true)
{
  cv::Mat imgFishAnterior,imgFishAnterior_bin,imgFishAnterior_Norm_tmplcrop;
  cv::Mat imgFishAnterior_Norm_bin_dense; //Used to Find Contours


 // Adapt Bounding Box size near image edges - Box Centrered at fishblob
  fishRotAnteriorBox.center = fishblob.pt;
  cv::Rect fishRotAnteriorBox_Bound = fishRotAnteriorBox.boundingRect();

  if (fishRotAnteriorBox_Bound.br().x > frame.cols) //Clip
        fishRotAnteriorBox.size.width = (int)(frame.cols - fishRotAnteriorBox.center.x)-1;
  if (fishRotAnteriorBox_Bound.tl().x < 0) // Shift
      fishRotAnteriorBox.center.x += abs(fishRotAnteriorBox_Bound.tl().x); //Shift Centre So BoundBox Fits in Frame

  //if ((fishblob.pt.y + fishRotAnteriorBox_Bound.size().height/2) > frame.rows)
  //    fishRotAnteriorBox.size.height = (int)(frame.rows - fishRotAnteriorBox.center.y)-1;
  if (fishRotAnteriorBox_Bound.br().y > frame.rows) //Clip
      fishRotAnteriorBox.center.y -= (int)(fishRotAnteriorBox_Bound.br().y-frame.rows)+1;

  if (fishRotAnteriorBox_Bound.tl().y < 0) // Shift
      fishRotAnteriorBox.center.y +=  abs(fishRotAnteriorBox_Bound.tl().y); //Shift Centre  so box To Fits in Frame


  fishRotAnteriorBox_Bound = fishRotAnteriorBox.boundingRect();

  if (fishRotAnteriorBox.size.width < gTrackerState.gszTemplateImg.width)
  {
        qDebug() << "scoreBlobRegion: fishRotAnteriorBox too narrow for templ";
        return -1;
  }
  if (fishRotAnteriorBox.size.height < gTrackerState.gszTemplateImg.height)
  {
     qDebug() << "scoreBlobRegion: fishRotAnteriorBox too short for templ";
     fishRotAnteriorBox.center.y = (int)(frame.rows -  gTrackerState.gszTemplateImg.height/2);
     fishRotAnteriorBox.size.height = gTrackerState.gszTemplateImg.height;
     return -1;
  }


  // Take bounding of blob,
  // get Rotated Box Centre Coords relative to the cut-out of the anterior Body - This we use to rotate the image
  // Allow Box Size to Adjust near Edges of Image
  //cv::RotatedRect fishRotAnteriorBox(fishblob.pt,
  //                                    cv::Size(std::min((int)fishblob.pt.x, gTrackerState.gFishBoundBoxSize),
  //                                             std::min((int)fishblob.pt.y,gTrackerState.gFishBoundBoxSize)),
  //                                              0);


   /// Size Of Norm Head Image

   cv::Size szFishAnteriorNorm = fishRotAnteriorBox_Bound.size();// (min(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height)+4,                              max(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height)+4);


   // Optimization to Make Search Region Adapt to size of Blob - Min is 2xtemplsize max is BoundBox Size
   //fishRotAnteriorBox_Bound.height = szFishAnteriorNorm.height = max(min(gTrackerState.gFishBoundBoxSize,szFishAnteriorNorm.height),2*gTrackerState.gszTemplateImg.height);
   //fishRotAnteriorBox_Bound.width = szFishAnteriorNorm.width = max(min(gTrackerState.gFishBoundBoxSize,szFishAnteriorNorm.width),2*gTrackerState.gszTemplateImg.height);
   //fishRotAnteriorBox_Bound.x = fishRotAnteriorBox.center.x - fishRotAnteriorBox_Bound.width/2; //Recenter Scaled Bound
   //fishRotAnteriorBox_Bound.y = fishRotAnteriorBox.center.y - fishRotAnteriorBox_Bound.height/2; //Recenter Scaled Bound

   // Define SCore Canvas - Where we draw results from Recognition Scanning
   cv::Mat maskRegionScore_Norm((int)szFishAnteriorNorm.height, (int)szFishAnteriorNorm.width, CV_32FC1, cv::Scalar(0.001, 0.001, 0.001));

   /// To Check Bounds Within Image
   cv::Rect imgBounds(0,0,frame.cols,frame.rows);

   // Check if region size is large enough to scan for fish
   if (!( //Looks Like a fish is found, now Check Bounds
       imgBounds.contains(fishRotAnteriorBox_Bound.br()) &&
           imgBounds.contains(fishRotAnteriorBox_Bound.tl())))
   {
            qDebug() << "scoreBlobRegion: boundingBox is out of frame bounds";
            return (-1); //This Fish Is out Of Bounds /
   }
   /// Check Minimum Size
   imgFishAnterior = frame(fishRotAnteriorBox_Bound); // getNormedBoundedImg(frame,fishRotAnteriorBox);


   if (szFishAnteriorNorm.area() < gTrackerState.szDNNClassifierImg.area())
   {
       //Not Enough Space Abort
       qDebug() << "scoreBlobRegion: Scan area too small for template size";
       return(0.0);
   }
  // Extract Region and rotate to orient larva body vertically
  // Use the FG Image to extract Head Frame
  //frame(fishRotAnteriorBox.boundingRect()).copyTo(imgFishAnterior);
  // Find point Center Coords of Image Region Of Search
  cv::Point2f ptRotCenter = fishblob.pt - fishRotAnteriorBox.boundingRect2f().tl();//fishRotAnteriorBox.center- fishRotAnteriorBox.boundingRect2f().tl();
  cv::Point ptmin;
  cv::Point ptmax = ptRotCenter; //Starting Assumption is that fish is located on Blob
  //Binarize Input To set Specific Sparseness/Density
  //imgFishAnterior_Norm_bin = sparseBinarize(imgFishAnterior_Norm,gTrackerState.fishnet_inputSparseness);
  imgFishAnterior.copyTo(imgFishAnterior_bin);
  //cv::imshow(std::string("BIN_N") + regTag,imgFishAnterior_Norm_bin);


  /// SliDing Window Scanning //- gTrackerState.gszTemplateImg.width/2
  /// \note for some reason X coord needs to be inverted so as to pick the right spot on imgFishAnterior_bin
  int iSlidePx_H_begin = std::max(0, (int)ptRotCenter.x - gTrackerState.szDNNClassifierImg.width/2-iSlidepxLim/2);// - iSlidepxLim/2 //max(0, imgFishAnterior_Norm.cols/2- sztemplate.width);
  int iSlidePx_H_lim = min(imgFishAnterior_bin.cols - gTrackerState.szDNNClassifierImg.width, iSlidePx_H_begin + iSlidepxLim);  //imgFishAnterior_Norm.cols/2; //min(imgFishAnterior_Norm.cols-sztemplate.width, max(0,imgFishAnterior_Norm.cols/2+ sztemplate.width) ) ;

   // V step - scanning for fishhead like image in steps
  //Bound Starting point
  int iSlidePx_V_begin = std::min(imgFishAnterior_bin.rows,
                                  std::max(0,(int)(ptRotCenter.y - gTrackerState.szDNNClassifierImg.height/2-iSlidepxLim/2) )); //  -  - iSlidepxLim/2==== -iSlidepxLim/2op(int)(ptRotCenter.y - sztemplate.height) sztemplate.height/2
  int iSlidePx_V_lim = min(imgFishAnterior_bin.rows-gTrackerState.szDNNClassifierImg.height, iSlidePx_V_begin + iSlidepxLim);//min(imgFishAnterior_Norm.rows - gTrackerState.gszTemplateImg.height, iSlidePx_V_begin + 10); //(int)(sztemplate.height/2)


  float scoreFish,scoreNonFish,scoreHuntMode,dscore,max_dscore = 0.0f; //Recognition Score tested in both Vertical Directions
  float mxLocFishScore,mxLocHuntScore = 0.0f; // Classification Scores for mode at point of best FishClass Match
  // Do netDetect using a Sliding window
  cv::Mat imgFishAnterior_Norm_tmplcrop_vflip;
  for (int i=iSlidePx_H_begin;i <= iSlidePx_H_lim;i+=iSlidePx_H_step)
  {
      for (int j=iSlidePx_V_begin;j <= iSlidePx_V_lim;j+=iSlidePx_V_step)
      {
          //Define Regions Top-Left Corner Position within Extracted Blob Region and
          cv::Point ptTopLeftTemplate(i,j);
          cv::Rect rectFishTemplateBound = cv::Rect(ptTopLeftTemplate, gTrackerState.szDNNClassifierImg  );
          //CROP Extract a Template sized subregion of Orthonormal Fish
          imgFishAnterior_bin(rectFishTemplateBound).copyTo(imgFishAnterior_Norm_tmplcrop);
          //dscore combines fish Score + huntMode score
          dscore = this->netDNNDetect_fish(imgFishAnterior_Norm_tmplcrop,scoreFish,scoreHuntMode,scoreNonFish);
          // Overide with Fish Score Only - As frequently Fish is missed
          dscore = max(scoreFish,scoreHuntMode); //Either Of these Being High
          ///Do Not Test Orientation - Blob Should Have the correct Angle
          //Check Both Vertical Orientations
          //cv::flip(imgFishAnterior_Norm_tmplcrop, imgFishAnterior_Norm_tmplcrop_vflip, 0);
          /// Score Number of Bin Pix In cropped region so as to push for regions that fully contain the eyes
          //float activePixRatio = (1+cv::sum(imgFishAnterior_Norm_tmplcrop)[0])/(imgFishAnterior_Norm_tmplcrop.cols*imgFishAnterior_Norm_tmplcrop.rows);
          ///  Store recognition score in Mask at(row,col) -//
          maskRegionScore_Norm.at<float>(j + gTrackerState.szDNNClassifierImg.height/2, // , row
                                         i + gTrackerState.szDNNClassifierImg.width/2) = dscore;//col //max(0.0f,dscore);//(scoreFish + scoreNonFish + 1e-3))/activePixRatio; //+ activePixRatio; //
          if (dscore > max_dscore && scoreNonFish < 0.6)
          {
              max_dscore = dscore;
              mxLocFishScore = scoreFish;
              mxLocHuntScore = scoreHuntMode;
              ptmax.y = j + gTrackerState.szDNNClassifierImg.height/2;
              ptmax.x = i + gTrackerState.szDNNClassifierImg.width/2;
          }

          if (bstopAtFirstMatch && max_dscore > gTrackerState.fishnet_classifier_thres)
          {
              //if (mxLocHuntScore > gTrackerState.fishnet_classifierHuntMode_thres)
                  //qDebug("Hunt Mode On");
               break;
          }
          //qDebug() << "(" << i+sztemplate.width/2 << "," <<j+sztemplate.height/2<<") = " << round(sc1*100)/100.0;
            //qDebug() << maskRegionScore_Norm.at<float>(j+gTrackerState.gszTemplateImg.height/2,
            //                                           i+gTrackerState.gszTemplateImg.width/2);
        }//For Each Vertical

      if (bstopAtFirstMatch && max_dscore > gTrackerState.fishnet_classifier_thres)
           break;
    }//For Each Horizontal

    /// Find and Best scorring point ScoreMask - Along with Normed Fish Image centered at best-point
    //Need to Rotate Score Image Back to Original Orientiation to return Coordinates oF Best Match
    //- Removed  - Scores scaled too low Find Max Match Position In Non-Norm pic (original orientation)
    double minL1, maxL1;
//    cv::GaussianBlur(maskRegionScore_Norm,maskRegionScore_Norm,
//                     cv::Size(max(2*(iSlidePx_H_step),2)+1,max(2*(iSlidePx_V_step),20)+1),
//                     max(iSlidePx_H_step/2,30),max(iSlidePx_V_step/2,30),BORDER_ISOLATED);

    cv::normalize(maskRegionScore_Norm, maskRegionScore_Norm, max_dscore, 0, cv::NORM_MINMAX);
    //ptmax = mxLocFishScore;
    cv::minMaxLoc(maskRegionScore_Norm,&minL1,&maxL1,&ptmin,&ptmax);
    //max_dscore = maxL1;
    // Rotate Max Point Back to Original Orientation
    //cv::Mat MrotInv = cv::getRotationMatrix2D( ptRotCenter, -fishblob.angle,1.0); //Rotate Upwarte Upwards
    //cv::warpAffine(maskRegionScore,outmaskRegionScore,MrotInv,szFishAnteriorNorm);
    //cv::Point ptmax_orig = rotateAboutPoint(ptmax,cv::Point2f(maskRegionScore_Norm.cols/2,maskRegionScore_Norm.rows/2),
    //                                        (fishblob.angle)*(CV_PI/180.0) ); //-fishblob.angle   angle
    // End Of FishNet Detection //


    //Update Blob Location And add Classifier Score if Higher than existing
     //Save Recognition Score - Don t use the Gaussian Blurred one -Too low
    if (max_dscore > fishblob.response)
    {
        fishblob.response = max_dscore;
        fishblob.FishClassScore = mxLocFishScore;
        fishblob.HuntModeClassScore = mxLocHuntScore;
        fishblob.pt = ptmax+fishRotAnteriorBox_Bound.tl(); //Shift Blob Position To Max  To Max Recognition Point
        //cv::circle(maskRegionScore_Norm,ptmax,3,CV_RGB(max_dscore,max_dscore,max_dscore),1);
    }else
    {
        //cv::circle(imgFishAnterior,ptmax,3,CV_RGB(250,250,250),1); //Indicate Max Score position By Grey Circle
        //cv::imshow("imgFishAnterior scoreBlobRegion "  + regTag, imgFishAnterior);
        //cv::imshow(("FishNet ScoreRegion (Norm) ") + regTag, maskRegionScore_Norm);
        //qDebug() << "scoreBlobRegion :" << max_dscore << " - Blob had higher class. score :" << fishblob.response ;
    }



    // DEBUG IMG //
    //cv::circle(imgFishAnterior,ptmax,4,CV_RGB(250,200,210),2);
    //cv::circle(frame,fishblob.pt,4,CV_RGB(250,200,210),2);

    //cv::imshow(string("Fish Region Body ") + regTag,imgFishAnterior);
    // DEBUG IMG //
    //cv::circle(imgFishAnterior_Norm,ptmax,3,CV_RGB(200,200,210),2);
    //cv::circle(maskRegionScore_Norm,ptmax,3,CV_RGB(0,0,0),2);


    //cv::imshow(string("Fish Region Body Norm ") + regTag,imgFishAnterior_Norm_bin);
    // DEBUG IMG //
    //cv::normalize(maskRegionScore_Norm, maskRegionScore_Norm, 0, 1, cv::NORM_MINMAX);
    //cv::imshow(string("Score Mask Body Norm") + regTag,maskRegionScore_Norm);                                   gionScore_Norm);

    /// Find Max Score Coords In Normed FishAnterior / Around Best Match Region (Using Normed Regiest Match Region (Using Normed Region)

    cv::Point ptTopLeftTemplate(max(0,ptmax.x-gTrackerState.szDNNClassifierImg.width/2-1),      //   TemplateImg.width/2),
                                max(0,ptmax.y-gTrackerState.szDNNClassifierImg.height/2-1) );
    // Stick To Boundary for Template Size Window
    ptTopLeftTemplate.x =((ptTopLeftTemplate.x + gTrackerState.szDNNClassifierImg.width) >= maskRegionScore_Norm.cols)?max(0,maskRegionScore_Norm.cols-gTrackerState.szDNNClassifierImg.width) : ptTopLeftTemplate.x;
    ptTopLeftTemplate.y =((ptTopLeftTemplate.y + gTrackerState.szDNNClassifierImg.height) >= maskRegionScore_Norm.rows)?max(0,maskRegionScore_Norm.rows-gTrackerState.szDNNClassifierImg.height) : ptTopLeftTemplate.y;
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


    cv::Rect rectFishTemplateBound = cv::Rect(ptTopLeftTemplate, gTrackerState.szDNNClassifierImg );
    assert(rectFishTemplateBound.width + rectFishTemplateBound.x <= imgFishAnterior.cols );
    assert(rectFishTemplateBound.height + rectFishTemplateBound.y <= imgFishAnterior.rows );

    /// CROP Extract a Template sized subregion of Orthonormal Fish ///
    imgFishAnterior_Norm_tmplcrop       = imgFishAnterior(rectFishTemplateBound);
    imgFishAnterior_Norm_tmplcrop.copyTo(outframeAnterior_Norm);
    maskRegionScore_Norm.copyTo(outmaskRegionScore); //(rectFishTemplateBound)

    /// Set Mark Point For Eye Detection ///
    //if (gTrackerState.bAdaptEyeMaskVOffset)
    //    gTrackerState.giHeadIsolationMaskVOffset = min(maxpt.y,imgFishAnterior_Norm.rows);//ptmax.y
    //    if (gTrackerState.bshowDetectorDebugImg)

    if (gTrackerState.bshowDetectorDebugImg ){
        cv::normalize(outmaskRegionScore, outmaskRegionScore, 1, 0, cv::NORM_MINMAX);
        cv::imshow(("FishNet ScoreRegion (Norm) ") + regTag, outmaskRegionScore);
        cv::imshow("imgFishAnterior scoreBlobRegion "  + regTag, imgFishAnterior);
    }



  return (max_dscore);

}


/// \brief Uses DNN classifier to detect most likely direction of fish.
/// \return best classifier output score/probability achieved while scanning rotations
/// @param regTag an Id for debugging purposes
/// @outframeAnterior_Norm returns image of isolated head centered at best detection point according to NN

float fishdetector::scoreBlobOrientation(cv::Mat frame,zftblob& fishblob,cv::Mat& outframeAnterior_Norm,
                                    cv::Mat& outmaskRegionScore,string regTag="0")
{
    cv::Mat imgFishAnterior_Norm,imgFishAnterior_Norm_best;
    float scoreFish,scoreNonFish,fRR,maxfRR = 0.0;

    int bestAngle = fishblob.angle;
    int startAngle = fishblob.angle;
    int dropStepCount =0; //Counts number consecutive scan points that the classifier signal drops So as to Stop Scanning early
    for (int a=(fishblob.angle-20);a<(startAngle+350);a+=2)
    {
        fishblob.angle = a;

        cv::RotatedRect fishRotAnteriorBox(fishblob.pt,  gTrackerState.gszTemplateImg,fishblob.angle);

         // To Check Bounds Within Image
         cv::Rect imgBounds(0,0,frame.cols,frame.rows);
         // Check if region size is large enough to scan for fish
         if (!( //Looks Like a fish is found, now Check Bounds
             imgBounds.contains(fishRotAnteriorBox.boundingRect().br()) &&
                 imgBounds.contains(fishRotAnteriorBox.boundingRect().tl())))
                  continue; //This Fish Is out Of Bounds /
        // Get Oriented image of fish Anterior
        imgFishAnterior_Norm =  getNormedTemplateImg(frame,fishRotAnteriorBox,false);
        // Check if Upright fish is detected within box
        fRR = netDNNDetect_normedfish(imgFishAnterior_Norm,scoreFish,scoreNonFish);

        if (fRR > maxfRR)
        {
            maxfRR = fRR;
            bestAngle = fishblob.angle;
            dropStepCount = 0; //Reset fail
            imgFishAnterior_Norm.copyTo( imgFishAnterior_Norm_best);
        }
        if (fRR < maxfRR) //If Moving away from peak/Down Gradient Stop Search
            dropStepCount++; //Break If Classifier threshold has been found

        if (dropStepCount > 10) //If Moving away from peak/Down Gradient Stop Search
            break;


    }// Test Full Circle

    fishblob.angle = (bestAngle)%360; //save best angle according to classifier (Convert from opencv Rotated Bound angle 0 being horizontal to tracker ref 0 on vertical

    if (maxfRR < gTrackerState.fishnet_classifier_thres) //No Match Found Across Angles
        fishblob.angle = (startAngle); //Reset - Angle Detection Failed
    else
    {
        //fishblob.response =  maxfRR;
        qDebug() << "Angle " << startAngle << "->" <<  fishblob.angle  << " fRR:" << maxfRR;
        cv::imshow("BestAngle",imgFishAnterior_Norm_best);
    }


return (maxfRR);

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

float fishdetector::netDNNDetect_fish(cv::Mat imgRegion_bin,float &fFishClass,float & fHuntModeClass,float & fNonFishClass)
{
    //std::cout << ".  run prediction..." << std::endl;
    cv::resize(imgRegion_bin,imgRegion_bin,{28,38});
    //cv::imshow("DNN detect",imgRegion_bin);
    ///\note if this crashes check input and output layer names match those of model
    std::vector< std::vector< float > > results = m_TFmodel_loc.predict<std::vector<float>>( {imgRegion_bin} );


    // print results - Assume Single Result Vector
    if ( results.size() > 0 )
    {
      //std::cout << "Output vector #" << i << ": ";
      for ( size_t j = 0; j < results[0].size(); j++ )
      {
        //qDebug() << results[0][j] << "\t";
         //std::cout << std::fixed << std::setprecision(4) << results[0][j] << "\t";
      }
       /// \TODO Exclude Small Fish From Here
       fFishClass = results[0][0]; //+results[0][2]; //Add Small and Large Fish Class Togetherresults[0][1] //Shifted
       fHuntModeClass = results[0][1];
       fNonFishClass = results[0][2];//1.0f - fFishClass;//results[0][2];
    }
    //std::cout << std::endl;
 //##Invert so Max Output predicts Fish
 //#return(1.0f-fFishClass);
   return((fFishClass+fHuntModeClass)-fNonFishClass);
}


///\deprecated An oriented classifier that attempts to detect Orientation of larva
float fishdetector::netDNNDetect_normedfish(cv::Mat imgRegion_bin,float &fFishClass,float & fNonFishClass)
{
    //std::cout << ".  run prediction..." << std::endl;
    cv::resize(imgRegion_bin,imgRegion_bin,{28,38});
    //cv::imshow("DNN detect",imgRegion_bin);
    std::vector< std::vector< float > > results = m_TFmodel_dir.predict<std::vector<float>>( {imgRegion_bin} );


    // print results - Assume Single Result Vector
    if ( results.size() > 0 )
    {
      //std::cout << "Output vector #" << i << ": ";
      for ( size_t j = 0; j < results[0].size(); j++ )
      {
        //qDebug() << results[0][j] << "\t";
         //std::cout << std::fixed << std::setprecision(4) << results[0][j] << "\t";
      }

       fFishClass = results[0][0];//+results[0][1]; //Fish + Fish SHifted
       fNonFishClass = 1.0f-fFishClass;// results[0][2];
    }
    //std::cout << std::endl;
//##Invert so Max Output predicts Fish
 return(1.0f-fFishClass);
}

///\brief Unit test fishnet detector by loading and testing against the set of images used during training (R Script)
///  and comparing output to the networks output to the script in R
void fishdetector::test()
{
    qDebug() << "Testing fishNet trained NN Classifier - loaded from YML Using stock template pics from disk.";
    QString strDirFish("/home/kostasl/workspace/zebrafishtrack/tensorDNN/test/fish/");
    QString strDirNonFish("/home/kostasl/workspace/zebrafishtrack/tensorDNN/test/nonfish/");

    std::vector<cv::Mat> vfish_mat = loadTemplatesFromDirectory(strDirFish);
    float fsumErrF =0.0f;
    float fsumErrNF =0.0f;
    float fCorrectF_class = 0.0;
    float fCorrectNF_class = 0.0;

    qDebug() << "~~~Test Fish templates~~~";
    float fishClassScore,nonfishScore,huntModeScore,dscore;
    for (int i=0;i<vfish_mat.size();i++)
    {
        cv::Mat imgTempl = vfish_mat[i];

        dscore = gTrackerState.fishnet.netDNNDetect_fish(imgTempl,fishClassScore,huntModeScore,nonfishScore);
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

        dscore = gTrackerState.fishnet.netDNNDetect_fish(imgTempl,fishClassScore,huntModeScore,nonfishScore);
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


