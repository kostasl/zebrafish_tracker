/// \brief This class will contain the learning and heurestic required to detect fish and manage the update of their respective model instances
///
///

//#include <config.h>
//#include "larvatrack.h"

//#include "fishmodel.h"
#include "fishdetector.h"
#include "template_detect.h"

extern trackerState gTrackerState;

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
    mW_L1 = loadImage(std::string("/home/kostasl/workspace/zebrafishtrack/Rplots/KC_SparseNet.pgm"));
    mW_L2 = loadImage(std::string("/home/kostasl/workspace/zebrafishtrack/Rplots/outputLayer_trained.pgm") );

    cv::threshold(mW_L1,mW_L1,0.1,1,cv::THRESH_BINARY);
    cv::threshold(mW_L2,mW_L2,0.1,1,cv::THRESH_BINARY);

    mW_L1.convertTo(mW_L1, CV_32FC1);
    mW_L2.convertTo(mW_L2, CV_32FC1);
}

/// \brief Uses Neural Net to scan region around provided blog and provide a detection score as a mask (returns max score value)
/// \todo could do image Pyramids to scan Across Scales
/// @param regTag an Id for debugging purposes
/// @outframeAnterior_Norm returns image of isolated head centered at best detection point according to NN
float fishdetector::scoreBlobRegion(cv::Mat frame,zftblob& fishblob,cv::Mat& outframeAnterior_Norm,cv::Mat& outmaskRegionScore,string regTag="0")
{
  cv::Mat imgFishAnterior,imgFishAnterior_Norm,imgFishAnterior_Norm_bin,imgFishAnterior_Norm_tmplcrop;
  // Take bounding of blob,
  //get Rotated Box Centre Coords relative to the cut-out of the anterior Body - This we use to rotate the image
   cv::RotatedRect fishRotAnteriorBox(fishblob.pt,
                                      cv::Size(2*gTrackerState.gFishBoundBoxSize,2*gTrackerState.gFishBoundBoxSize),
                                                fishblob.angle);
   /// Size Of Norm Head Image
   cv::Size szFishAnteriorNorm = fishRotAnteriorBox.boundingRect().size();// (min(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height)+4,                              max(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height)+4);
   // Define SCore Canvas - Where we draw results from Recognition Scanning
   cv::Mat maskRegionScore_Norm((int)szFishAnteriorNorm.width, (int)szFishAnteriorNorm.height, CV_32FC1, cv::Scalar(0, 0, 0));

   //To Check Bounds Within Image
   cv::Rect imgBounds(0,0,frame.cols,frame.rows);
   cv::Size sztemplate = gTrackerState.gLastfishimg_template.size();

   // Check if region size is large enough to scan for fish
   if (!( //Looks Like a fish is found, now Check Bounds
       imgBounds.contains(fishRotAnteriorBox.boundingRect().br()) &&
           imgBounds.contains(fishRotAnteriorBox.boundingRect().tl())))
            return (-1); //This Fish Is out Of Bounds /

  /// Extract Region and rotate to orient larva body vertically
  // Use the FG Image to extract Head Frame
  frame(fishRotAnteriorBox.boundingRect()).copyTo(imgFishAnterior);

  /// Make Rotation MAtrix About Centre Of Cropped Image
  cv::Point2f ptRotCenter = fishRotAnteriorBox.center - fishRotAnteriorBox.boundingRect2f().tl();
  cv::Mat Mrot = cv::getRotationMatrix2D( ptRotCenter, fishblob.angle,1.0); //Rotate Upwards
  ///Make Rotation Transformation
  //Need to fix size of Upright/Normed Image
  cv::warpAffine(imgFishAnterior,imgFishAnterior_Norm,Mrot,szFishAnteriorNorm);

  // Binarize Through Adaptive Threshold to enhance fish-like pattern // Substract - const val from mean
  // Ideally We want to maintain input sparseness
  double activePixRatio = 1.0;
  int meanAdapt = 0;
  int maxIter = 20;
  // Regulate Input Sparseness
  while (activePixRatio > gTrackerState.fishnet_inputSparseness &&
         maxIter > 0)
  {
    cv::adaptiveThreshold(imgFishAnterior_Norm,imgFishAnterior_Norm_bin,1,cv::ADAPTIVE_THRESH_MEAN_C,cv::THRESH_BINARY,szFishAnteriorNorm.width,meanAdapt);
    activePixRatio = (cv::sum(imgFishAnterior_Norm_bin)[0])/(imgFishAnterior_Norm_bin.cols*imgFishAnterior_Norm_bin.rows);
    meanAdapt += (int)255*(gTrackerState.fishnet_inputSparseness-activePixRatio)-2;
    maxIter--;
  }
  //cv::threshold(imgFishAnterior_Norm,imgFishAnterior_Norm_bin,gTrackerState.g_Segthresh,1,cv::THRESH_BINARY);
  if (!imgFishAnterior_Norm_bin.empty())
      cv::imshow("FIshBody Norm Bin"+regTag,imgFishAnterior_Norm_bin*255);

  //cv::Point ptTopLeftTemplate(szFishAnteriorNorm.width/2-gTrackerState.gLastfishimg_template.size().width/2,
  //                         szFishAnteriorNorm.height/2-gTrackerState.gLastfishimg_template.size().height/2);

   //Check Center Of Image / If Something is found Do Sliding Window

  /// SliDing Window Scanning
  int iSlidePx_H_begin = ptRotCenter.x- sztemplate.width/2-5;//max(0, imgFishAnterior_Norm.cols/2- sztemplate.width);
  int iSlidePx_H_lim = iSlidePx_H_begin+10;  //imgFishAnterior_Norm.cols/2; //min(imgFishAnterior_Norm.cols-sztemplate.width, max(0,imgFishAnterior_Norm.cols/2+ sztemplate.width) ) ;
    int iSlidePx_H_step = 3;

  int iSlidePx_V_begin = std::max(0,(int)(ptRotCenter.y - sztemplate.height/2)-5); //(int)(ptRotCenter.y - sztemplate.height) sztemplate.height/2
  int iSlidePx_V_lim = min(imgFishAnterior_Norm.rows - sztemplate.height, iSlidePx_V_begin +10); //(int)(sztemplate.height/2)
  int iSlidePx_V_step = 3;

  float scoreFish,scoreNonFish,dscore; //Recognition Score tested in both Vertical Directions
  // Do netDetect using a Sliding window
  cv::Mat imgFishAnterior_Norm_tmplcrop_vflip;
  for (int i=iSlidePx_H_begin;i <= iSlidePx_H_lim;i+=iSlidePx_H_step)
  {
      for (int j=iSlidePx_V_begin;j <= iSlidePx_V_lim;j+=iSlidePx_V_step)
      {
          //Define Regions Top-Left Corner Position within Extracted Blob Region and
          cv::Point ptTopLeftTemplate(i,j);
          cv::Rect rectFishTemplateBound = cv::Rect(ptTopLeftTemplate, sztemplate );
          //CROP Extract a Template sized subregion of Orthonormal Fish
          imgFishAnterior_Norm_bin(rectFishTemplateBound).copyTo(imgFishAnterior_Norm_tmplcrop);

          dscore = this->netDetect(imgFishAnterior_Norm_tmplcrop,scoreFish,scoreNonFish);
          ///Do Not Test Orientation - Blob Should Have the correct Angle
          //Check Both Vertical Orientations
          //cv::flip(imgFishAnterior_Norm_tmplcrop, imgFishAnterior_Norm_tmplcrop_vflip, 0);

          /// Score Number of Bin Pix In cropped region so as to push for regions that fully contain the eyes
          float activePixRatio = (1+cv::sum(imgFishAnterior_Norm_tmplcrop)[0])/(imgFishAnterior_Norm_tmplcrop.cols*imgFishAnterior_Norm_tmplcrop.rows);
          ///  Store recognition score in Mask at(row,col) -//
          maskRegionScore_Norm.at<float>(j+sztemplate.height/2, i+sztemplate.width/2) = (scoreFish/(scoreFish + scoreNonFish + 1e-3))/activePixRatio; //+ activePixRatio; //

          //qDebug() << "(" << i+sztemplate.width/2 << "," <<j+sztemplate.height/2<<") = " << round(sc1*100)/100.0;

        }//For Each Vertical
    }//For Each Horizontal

    /// Find and Best scorring point ScoreMask - Along with Normed Fish Image centered at best-point
    //Need to Rotate Score Image Back to Original Orientiation to return Coordinates oF Best Match
    /// Find Max Match Position In Non-Norm pic (original orientation)
    double minL1,maxL1;
    cv::Point ptmin,ptmax;
    cv::GaussianBlur(maskRegionScore_Norm,maskRegionScore_Norm,cv::Size(9,9),15,15);
    cv::minMaxLoc(maskRegionScore_Norm,&minL1,&maxL1,&ptmin,&ptmax);

    // Rotate Max Point Back to Original Orientation
    //cv::Mat MrotInv = cv::getRotationMatrix2D( ptRotCenter, -fishblob.angle,1.0); //Rotate Upwards
    //cv::warpAffine(maskRegionScore,outmaskRegionScore,MrotInv,szFishAnteriorNorm);
    cv::Point ptmax_orig = rotateAboutPoint(ptmax,cv::Point2f(maskRegionScore_Norm.cols/2,maskRegionScore_Norm.rows/2),
                                            (fishblob.angle)*(CV_PI/180.0) ); //-fishblob.angle

    //Update Blob Location And add Classifier Score
    fishblob.response = maxL1; //Save Recognition Score
    fishblob.pt = ptmax_orig+fishRotAnteriorBox.boundingRect().tl(); //Shift Blob Position To Max Recognition Point

    // DEBUG IMG //
    cv::circle(imgFishAnterior,ptmax_orig,3,CV_RGB(250,200,210),2);
    //cv::imshow(string("Fish Region Body ") + regTag,imgFishAnterior);
    // DEBUG IMG //
    cv::circle(imgFishAnterior_Norm,ptmax,3,CV_RGB(250,200,210),2);
    //cv::imshow(string("Fish Region Body Norm ") + regTag,imgFishAnterior_Norm);
    // DEBUG IMG //
    cv::normalize(maskRegionScore_Norm, maskRegionScore_Norm, 0, 1, cv::NORM_MINMAX);
    //cv::imshow(string("Score Mask Body Norm") + regTag,maskRegionScore_Norm);

    /// Find Max Score Coords In Normed FishAnterior / Around Best Match Region (Using Normed Region)
    cv::Point ptTopLeftTemplate(szFishAnteriorNorm.width/2-gTrackerState.gLastfishimg_template.size().width/2,
                             szFishAnteriorNorm.height/2-gTrackerState.gLastfishimg_template.size().height/2);

//    // Bound search region so cropping Of Norm does not exceed image bounds
//    cv::Rect rect_bound = cv::Rect(cv::Point(gTrackerState.gLastfishimg_template.size().width/2,
//                                             gTrackerState.gLastfishimg_template.size().height/2),
//                                             cv::Size(maskRegionScore.cols-gTrackerState.gLastfishimg_template.size().width/2,
//                                                       maskRegionScore.rows-gTrackerState.gLastfishimg_template.size().height/2));
//    cv::minMaxLoc(maskRegionScore(rect_bound),&minL1,&maxL1,&ptmin,&ptmax);
//    cv::Point ptTopLeftTemplate(min(imgFishAnterior_Norm.cols, max(0, ptmax.x - gTrackerState.gLastfishimg_template.size().width/2)),
//                                min(imgFishAnterior_Norm.rows, max(0, ptmax.y - gTrackerState.gLastfishimg_template.size().height/2))
//                               );


    cv::Rect rectFishTemplateBound = cv::Rect(ptTopLeftTemplate, sztemplate );
    assert(rectFishTemplateBound.width + rectFishTemplateBound.x <= imgFishAnterior_Norm.cols );
    assert(rectFishTemplateBound.height + rectFishTemplateBound.y <= imgFishAnterior_Norm.rows );

    /// CROP Extract a Template sized subregion of Orthonormal Fish ///
    imgFishAnterior_Norm_tmplcrop       = imgFishAnterior_Norm(rectFishTemplateBound);
    imgFishAnterior_Norm_tmplcrop.copyTo(outframeAnterior_Norm);
    maskRegionScore_Norm(rectFishTemplateBound).copyTo(outmaskRegionScore);

  return (maxL1);

}

/// \brief Applies pre-trained MB like NN on Binarized Input image
/// Networks supports two L2 neurons - These recognition nets suffer from decreases in input sparseness :
/// More active inputs increase the output scores - reducing the networks selectivity
float fishdetector::netDetect(cv::Mat imgRegion_bin,float &fFishClass,float & fNonFishClass)
{
    fL1_activity_thres = gTrackerState.fishnet_L1_threshold;

    // Input Is converted to Row Vector So we can do Matrix Multiplation
    assert(imgRegion_bin.cols*imgRegion_bin.rows == mW_L1.rows);
    cv::Mat vIn = imgRegion_bin.reshape(1,mW_L1.rows).t();  //Row Vector
    vIn.convertTo(vIn, CV_32FC1);

    // operation multiplies matrix A of size [a x b] with matrix B of size [b x c]
    //to the Layer 1 output produce matrix C of size [a x c]
    mL1_out = vIn*mW_L1;

    // Threshold for Activation Function
    cv::threshold(mL1_out,mL1_out,fL1_activity_thres,1,cv::THRESH_BINARY);
    //Calc Output
    mL2_out =  mL1_out*mW_L2;

    //cv::imshow("L1 Out", mL1_out);
    //Output fraction of Active Input that is filtered by Synaptic Weights, (Fraction of Active Pass-through KC neurons)
    fFishClass = mL2_out.at<float>(0,0)/mW_L1.cols;
    // Check 2 row (neuron) output
    fNonFishClass = mL2_out.at<float>(0,1)/mW_L1.cols;

    //double minL1,maxL1;
    //cv::minMaxLoc(mL1_out,&minL1,&maxL1);
    //qDebug() << "***R: " << fOut << " KCmin: "<< minL1 << " KCmax: " << maxL1;
    //cv::imshow("output",mL2_out);

    return(fFishClass-fNonFishClass);
}

