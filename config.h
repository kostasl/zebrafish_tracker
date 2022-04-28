#ifndef CONFIG_H
#define CONFIG_H


#define ZTF_FISHCONTOURSIZE          62//40
#define ZTF_TAILFITMAXITERATIONS     200 //For Spine To Contour Tail Fitting
#define ZTF_TAILSPINECOUNT          8
#define EYE_SEG_SAMPLE_POINTS_COUNT 20

#if defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    #define USE_CUDA
#else
    #undef USE_CUDA
#endif

#undef  USE_CUDA_FOR_TEMPLATE_MATCHING //Can Be Very Slow When Using Current Algorithm for searching through Cache
#define TEMPLATE_COL_SEARCH_REGION  15 //Optimization How Far to Search Along Col in the Cache From the Starting POint

#define G_MAX_FOOD_BLOB_LIMIT = 150


#include <limits>
#include <string>
#include <cmath>
#include <QString>
#include <string.h>

#include <errorhandlers.h> // My Custom Mem Fault Handling Functions and Debug
//#include "eyesdetector.h"
#include <ltROI.h> //Defines the ROI types
//#include <fishmodel.h>
//#include <foodmodel.h>
#include <zfttracks.h>
#include <CSS/CurveCSS.h>
#include <template_detect.h>
#include <fishdetector.h>

//QT
#include <QDirIterator>
#include <QDir>
#include <QDebug>
//#include <QThread>
#include <QTime>
#include <QFileDialog>
#include <QElapsedTimer>

//Open CV
#include <opencv2/opencv_modules.hpp> //THe Cuda Defines are in here
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
//#include <opencv2/bgsegm.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/video/background_segm.hpp>

#include <opencv2/core/ocl.hpp> //For setting setUseOpenCL

/// For State Saving //
#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/common.hpp>
#include <fstream>

//Random
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/// Fish Detection CV BG Model //
extern cv::Ptr<cv::BackgroundSubtractorMOG2> pMOG2; //MOG2 Background subtractor


/// Contour Shaping - Global Variables
extern double sigma ;
extern int M ; //Gaussian Kernel Size

/// Morphological Kernels
extern cv::Mat kernelOpen;
extern cv::Mat kernelDilateMOGMask;
extern cv::Mat kernelOpenfish;
extern cv::Mat kernelClose;

extern cv::Mat gFishTemplateCache; //A mosaic image contaning copies of template across different angles


/// \brief pointPairs is a vector holding points that the user chooses so as to measure distances between arbitrary points
typedef std::pair<cv::Point,cv::Point>  pointPair;
typedef std::vector<pointPair > pointPairs;

// Gaussian Curve Smoothing Kernels For fish Contour//
//extern std::vector<double> gGaussian,dgGaussian,d2gGaussian;

extern QElapsedTimer gTimer;
extern QFile outfishdatafile;
extern QFile outfooddatafile;
extern QFile EyeDetectorRL; //Reinforcement Learned Behaviour For Eye Segmentation -

extern std::ofstream foutLog;//Used for Logging To File
static jmp_buf env; // For Skipping SIG sev and jumping to Prog Counter

/// Global Matrices Used to show debug images
extern cv::Mat frameDebugA,frameDebugB,frameDebugC,frameDebugD;
extern cv::Mat gframeCurrent,gframeLast; //Updated in processVideo Global Var Holding Copy of current and previous frame - usefull for opticflows
extern cv::Mat gframeBGImage;

/// \todo using a global var is a quick hack to transfer info from blob/Mask processing to fishmodel / Need to change the Blob Struct to do this properly
extern cv::Point gptHead,gptTail; //Candidate Fish Contour Position Of HEad - Use for template Detect

//extern ltROIlist vRoi;
//Rect Roi Keep Away from L-R Edges to Avoid Tracking IR lightRing Edges
//extern  cv::Point ptROI1 ;
//extern  cv::Point ptROI2;
//extern cv::Point ptROI3 ;
//extern cv::Point ptROI4 ;


class trackerState
{
 public:
     trackerState();
     void saveState(std::string strFilename);
     void loadState(std::string strFilename);
     void setVidFps(float fps);
     /// \brief Process user provided config params and set corresponding internal/global variables
     void initGlobalParams(cv::CommandLineParser& parser,QStringList& inVidFileNames); //Read Command Line/Config Options
     /// \brief Initializes ROI at start of tracking depending on user params / either large circle or user defined/configurable polygon
     void  initROI(uint framewidth,uint frameheight);

     /// \brief Load Q Resources
     static void loadFromQrc(QString qrc,cv::Mat& imRes,int flag = cv::IMREAD_COLOR); //Load Resources

      enum state {PAUSED,TRACKING,DIST_MEASURE,SAVING,EXITING};

      /// VIDEO AND BACKGROUND PROCESSING //
      float gfVidfps                   = 400;
      uint frame_pxwidth               = 640; //Video Frame pixel Dimensions/ Default Changed when Video Is opened
      uint frame_pxheight              = 480;

      const unsigned int MOGhistory   = 15; //Use 100 frames Distributed across the video length To Find What the BGModel is
      double gdMOGBGRatio             = 0.05; ///If a foreground pixel keeps semi-constant value for about backgroundRatio*history frames, it's considered background and added to the model as a center of a new component.
      double dBGMaskAccumulateSpeed             = 1.0/(4.0*MOGhistory);

      const double dLearningRate                = 1.0/(MOGhistory); //Learning Rate During Initial BG Modelling - Learn Slow So 1st Playbacl Frame doesnt look new anymore
      const double dLearningRateNominal         = 0.001; //Fast Rate as BG Learning Allows for threshold after BGSubstract operation to Augment The mask
      double dactiveMOGLearningRate             = dLearningRateNominal;
      //Processing Loop delay
      uint cFrameDelayms                        = 1;

      /// Constants ///
      /// BLOB DETECTION Filters //
      //Area Filters
      const double dMeanBlobArea                  = 100; //Initial Value that will get updated
      const double dVarBlobArea                   = 20;
      const unsigned int gc_fishLength            = 100; //px Length Of Fish
      const unsigned int thresh_fishblobarea      = 400; //Min area above which to Filter The fish blobs
      const unsigned int thresh_maxfishblobarea     = 3650; //Min area above which to Filter The fish blobs
      const unsigned int gthres_maxfoodblobarea     = thresh_fishblobarea/3;

      const int gFitTailIntensityScanAngleDeg   = 60; //
      const int gFishTailSpineSegmentCount      = ZTF_TAILSPINECOUNT;
      const int gcFishContourSize               = ZTF_FISHCONTOURSIZE;
      const int gMaxFitIterations               = ZTF_TAILFITMAXITERATIONS; //Constant For Max Iteration to Fit Tail Spine to Fish Contour

      int gcMaxFishModelInactiveFrames          = gfVidfps/2; //Number of frames inactive until track is deleted
      int gcMaxFoodModelInactiveFrames          = gfVidfps/5; //Number of frames inactive (Not Matched to a Blob) until track is deleted
      int gcMinFoodModelActiveFrames            = gfVidfps/20; //Min Number of consecutive frames it needs to be active  otherwise its deleted
      float gMaxClusterRadiusFoodToBlob           = 3; //Per Sec / This changes depending on FPS (setFPS)
      const int thActive                            = 0;// Deprecated If a track becomes inactive but it has been active less than thActive frames, the track will be deleted.
      const int gc_FishTailSpineSegmentLength_init  = 16;
      int gFoodReportInterval                       = (int)gfVidfps;
      const int nTemplatesToLoad                    = 11; //Number of Templates To Load Into Cache - These need to exist as images in QtResources
      const int gi_FoodModelNumberLimit             = 250; // Maximum Number of Food Objects /Prey To track

      int keyboard; //input from keyboard
      int screenx,screeny;

      cv::Mat gLastfishimg_template;
      //Rect/Polygon Roi Keep Away from L-R Edges to Avoid Tracking IR lightRing Edges
      std::shared_ptr<std::vector<tROIpt>> userROI;

      std::string  gstroutDirCSV;//Tracker's Output Directory
      std::string  gstroutDirTemplates;//Tracker's Source Of Larva Head Image Templates Library
      std::string  gstrinDirVid;
      QStringList inVidFileNames; //List of Video Files to Process
      std::string  gstrvidFilename; //Currently Tracked Vid
      cv::Size gszTemplateImg = cv::Size(28,38);
      cv::Size sztemplateArray_Icon=  cv::Size(std::max(gszTemplateImg.width, gszTemplateImg.height),
                                     std::max(gszTemplateImg.width, gszTemplateImg.height)
                                    );


      cv::Rect rect_pasteregion;//For Inset Pasting

      //Contours Shaping
      double curveSmoothKernelSigma;
      int curveSmoothKernelSize_M; //Gaussian Kernel Size

      uint gi_MaxFishID; // Records the state of Fish ID so as to increment unique ids
      uint gi_MaxFoodID; //Declared in Model Header Files

      QString outfilename;
      std::string gstrwinName = "Fishtracker Frame";

     // Global Control Vars ///
     /// \brief bTracking
     ///// Option Flags //
      bool bAllowOnlyOneTrackedItem = false;
      bool bshowMask                = false; //Debug option True will show the BGSubstracted IMage/Processed Mask
      bool bshowDetectorDebugImg    = false; //Debug option  True will show the classifier scoring Masks and Extracted Fish Anterior Images

      bool bStartPaused             = false; //Command line controlled
      bool bPaused                  = false;
      bool bExiting;
      bool bTracking          = true;
      bool bTrackFood         = true;
      //bool bSaveBlobsToFile   = false; //Check in fnct processBlobs - saves output CSV

      bool bAddPreyManually   = false;
      bool bMeasure2pDistance = true; /// A mode allowing 2point distance measurement
      bool bTrackFish         = true;
      bool bOnlyTrackFishinROI = false; /// Whether to attempt tracking outside of Defined ROI - or Stop Tracking of Fish when OUtside ROI,
      bool bTrackAllPrey      = false; ///Track All detected Prey Models/Once established Active
      bool bRecordToFile      = true;
      bool bSaveImages            = false;
      bool gOptimizeShapeMatching = false; ///Set to false To disable matchShapes in FindMatching Contour

      //GUI User Interactions
      bool b1stPointSet            = false;
      bool bMouseLButtonDown       = false;
      bool bROIChanged             = false;

      bool bEyesDetected                  = false; ///Flip True to save eye shape feature for future detection
      bool bStoreThisTemplate             = false;
      bool bDraggingTemplateCentre        = false;

      bool bStartFrameChanged         = false; /// When True, the Video Processing loop stops /and reloads video starting from new Start Position

      bool bRenderToDisplay           = true; ///Updates Screen to User When True
      bool bRenderWithAlpha           = false;
      bool bDrawFoodBlob              = true; ///Draw circle around identified food blobs (prior to model matching)
      bool bOffLineTracking           = false; ///Skip Frequent Display Updates So as to  Speed Up Tracking
      bool bBlindSourceTracking       = false; /// Used for Data Labelling, so as to hide the data source/group/condition
      bool bStaticBGMaskRemove        = false; /// Remove Pixs from FG mask that have been shown static in the Accumulated Mask after the BGLearning Phase
      bool bUseBGModelling                      = true; ///Use BG Modelling TO Segment FG Objects
      bool gbUpdateBGModel                      = true; //When Set a new BGModel Is learned at the beginning of the next video
      bool gbUpdateBGModelOnAllVids             = true; //When Set a new BGModel Is learned at the beginning of the next video
      bool bApplyFishMaskBeforeFeatureDetection = true; ///Pass the masked image of the fish to the feature detector /Fails If the Mask draw contour only has the edges
      bool bUseTemplateMatching                 = true; /// Cmd Line Param Use Template Matching Following DNN classifier Success
      bool bFitSpineToTail                      = true; // Periodically Runs The Contour And Tail Fitting Spine Optimization Algorith
      bool bUseContourToFitSpine                = true; // Periodically Runs The Contour And Tail Fitting Spine Optimization Algorith
      bool bSkipExisting                        = false; /// If A Tracker DataFile Exists Then Skip This Video
      bool bMakeCustomROIRegion                 = false; /// Uses Point array to construct
      bool bUseMaskedFishForSpineDetect         = true; /// When True, The Spine Is fit to the Masked FG Fish Image and not the full frame- (Masks can lose fine features)
      bool bTemplateSearchThroughRows           = true; /// Stops TemplateFind to Scan Through All Rows (diff template images)- speeding up search + fail - Rows still Randomly Switch between attempts
      bool bRemovePixelNoise                    = false; //Run Gaussian Filter Noise Reduction During Tracking
      bool bUseGPU                              = false;
      bool bUseOpenCL                           = true;
      bool bUseHistEqualization                 = true; //To enhance to contrast in Eye Ellipse detection
      bool bUseEllipseEdgeFittingMethod         = true; //Allow to Use the 2nd Efficient Method of Ellipsoid Fitting if the 1st one fails - Set to false to Make trakcing Faster
      bool bAdaptEyeMaskVOffset                 = true; // Check in fishDetector.cpp

      /// \todo Make this path relative or embed resource
      //string strTemplateImg = "/home/kostasl/workspace/cam_preycapture/src/zebraprey_track/img/fishbody_tmp.pgm";
      std::string strTemplateImg = ":/img/fish_template_"; ///Load From Resource

      ///Specific To the Tracked Video Options//
      uint uiCurrentFrame = 1;
      uint uiStartFrame = 1;
      uint uiStopFrame = 0;
      uint iSpineContourFitFramePeriod         = 20; //Check that Tail Fitting Matches Contour Every X Frames

      /// Segmentation / threshold  Params
      int g_FGSegthresh;            //Default Set By Command Line Params
      int g_FGStaticMaskSegthresh   = g_FGSegthresh;//Applied On THe BG Substracted Image / Image Threshold to segment BG - Fish Segmentation uses a higher 2Xg_Segthresh threshold
      int g_SegFoodThesMin        = max(5,g_FGSegthresh-25); //Low thres For Food Detection / Doing Gradual Step Wise with SimpleBlob
      int g_SegFoodThesMax        = g_SegFoodThesMin+15; //Up thres Scan For Food Detection / Doing Gradual Step Wise with SimpleBlob
      //int g_SegInnerthreshMult    = 3; //Image Threshold for Inner FIsh Features //Deprecated
      //int g_BGthresh              = 5; //BG threshold segmentation
      int gi_ThresholdMatching    = 10; /// Minimum Score to accept that a contour has been found

      /// Eye Tracking Params
      int gi_CannyThres           = 150;
      int gi_CannyThresSmall      = 50; //Aperture size should be odd between 3 and 7 in function Canny
      int gi_maxEllipseMajor      = 30; /// thres  for Eye Ellipse Detection methods
      int gi_minEllipseMajor      = 13; ///thres for Eye Ellipse Detection methods (These Values Tested Woodrked Best)
      int gi_minEllipseMinor      = 0; /// ellipse detection width - When 0 it allows for detecting straight line
      int gi_MaxEllipseSamples    = 10; //The number of fitted ellipsoids draw from the ranked queue to calculate mean fitted eye Ellipse
      int gi_VotesEllipseThres            = 5; //Votes thres for The Backup Ellipse Detection Based on the Hough Transform
      int gthresEyeSeg                    = -30; //-23 Additional Adjustment for Adaptive Threshold  For Eye Segmentation In Isolated Head IMage -Shown On GUI

      int gthresEyeSegL                   = 2;
      int gFishTailSpineSegmentLength     = 9;
      // Eye Masks //
      int iEyeHMaskSepRadius                = 18; //Radius of Mask centred at bottom of head, also used as Threshold Sampling Arc in Fish Head Mask
      //int giEyeIsolationMaskRadius       = 17; Not Used //Mask circle between eyes
      int iEyeVMaskSepWidth               = 25; //5 px width vertical line separates the eyes for segmentation
      int iEyeVMaskSepHeight              = 46;
      int eyeMaskVLineThickness           = 8; //Width Vertical Midline Separating The eyes

      /// Fishnet Classifier params //
      //float fishnet_L1_threshold  = 0.5; //L1 neuron Activity Threshold Sets the Pattern Selectivity and sparseness of L1 output
      float fishnet_classifier_thres  = 0.94f; //L1 neuron Activity Threshold Sets the Pattern Selectivity and sparseness of L1 output
      float fishnet_inputSparseness = 0.1f; //Ratio of Active Pixels in Binarized input Image

      // BackProp YAML model - DEpecrated
      /// \todo Make this Relative Paths - Adjustable
      string strBackPropModelYMLFile        = "/home/meyerlab/workspace/zebrafishtrack/Rplots/fishNet.yml";
      string strDNNTensorFlowModelFile      = "/home/meyerlab/workspace/zebrafishtrack/tensorDNN/savedmodels/fishNet_loc/";
      std::string strDNNTensorFlowVerticalModelFile = "/home/meyerlab/workspace/zebrafishtrack/tensorDNN/savedmodels/fishNet_dir/";

      ///Fish Features Detection Params
      int gFishTemplateAngleSteps     = 1;
      int gEyeTemplateAngleSteps      = 5;

      double eyeStepIncrement         = 0.8; //Eye Angles Can be Slowly Updated on each Frame- Change with Step Size eyeStepIncrement
      double gTemplateMatchThreshold  = 0.73; //Template Matching is tested After Fish Net Classifier Has passed-
      double gTemplateMatchThreshold_LowLimit = 0.65;
      double gTemplateMatchThreshold_UpLimit = 0.95;

      int gFishBoundBoxSize               = 100; ///100 For HRes Top CamB 24/ pixel width/radius of bounding Box When Isolating the fish's head From the image
      int gnumberOfTemplatesInCache       = 0; //INcreases As new Are Added
      float  gDisplacementThreshold       = 2.0; //Distance That Fish Is displaced so as to consider active and Record A point For the rendered Track /
      int  gDisplacementLimitPerFrame    = gFishBoundBoxSize*4; //Distance That Fish-Blob can be allowed to displace - Filter Out Large Motion Noise in FishModel UpdateState
      int  gAngleChangeLimitPerFrame    = 90; //Distance That Fish-Blob can be allowed to displace - Filter Out Large Motion Noise in FishModel UpdateState

      int iLastKnownGoodTemplateRow   = 0;
      int iFishAngleOffset            = 0;
      double gUserReward              = 0; //User feedback for reinforcement learning
      int iTemplateMatchFailCounter   = 0; //Counts the number of consecutive times template failed to match
      //using namespace std;

      /// Random Number Generator
      const gsl_rng_type * T;
      gsl_rng * p_gsl_r;

      // Other fonts:
      //   CV_FONT_HERSHEY_SIMPLEX, CV_FONT_HERSHEY_PLAIN,
      //   CV_FONT_HERSHEY_DUPLEX, CV_FONT_HERSHEY_COMPLEX,
      //   CV_FONT_HERSHEY_TRIPLEX, CV_FONT_HERSHEY_COMPLEX_SMALL,
      //   CV_FONT_HERSHEY_SCRIPT_SIMPLEX, CV_FONT_HERSHEY_SCRIPT_COMPLEX
      const int trackFnt = cv::FONT_HERSHEY_SIMPLEX;  //Font for Reporting - Tracking
      const float trackFntScale = 0.4f;
      /// Contour Shaping Gaussian Kernels //
      std::vector<double> gGaussian,dgGaussian,d2gGaussian;
      double dGaussContourKernelSigma = 6.0;
      int dGaussContourKernelSize = static_cast<int>(round((18.0*dGaussContourKernelSigma +1.0) / 2.0)) * 2 - 1; //Gaussian Kernel Size

      // List of ROIs
      ltROIlist vRoi;
      cv::Point ptROI1,ptROI2,ptROI3,ptROI4;
      int iROIRadius               = 320;

      //list of template images
      std::vector<cv::Mat> vTemplImg;

      fishdetector fishnet;

      cv::Mat mMOGMask;
      cv::Mat mfgFishMask;
      cv::Mat mfgFishFrame;
      cv::Mat mfgPreyMask;
    /// overload size operator / return full state object size
     size_t size() const _GLIBCXX_NOEXCEPT
    { return 1; //mStateValue.size()*mStateValue[1].size()*mStateValue[1][1].size();
     }

 // This method lets cereal know which data members to serialize
   template<class Archive>
   void serialize(Archive & archive)
   {
     archive(gstrwinName,CEREAL_NVP(gstroutDirCSV),CEREAL_NVP(gstrinDirVid),gstrvidFilename,outfilename.toStdString(),
            userROI,bRecordToFile,bTrackFish,bSaveImages,bUseEllipseEdgeFittingMethod,
            bTemplateSearchThroughRows,bApplyFishMaskBeforeFeatureDetection,bUseOpenCL,bUseGPU,bBlindSourceTracking,bStaticBGMaskRemove,
            gbUpdateBGModel,gbUpdateBGModelOnAllVids,bFitSpineToTail,bUseMaskedFishForSpineDetect,bUseHistEqualization,bRemovePixelNoise,bMeasure2pDistance,bAddPreyManually,
            bRenderToDisplay,bRenderWithAlpha,bOffLineTracking,bDrawFoodBlob,curveSmoothKernelSigma,curveSmoothKernelSize_M,
            bUseBGModelling,gbUpdateBGModel,gbUpdateBGModelOnAllVids,
            bSkipExisting,bTrackFood,bTracking,bStartPaused,uiStartFrame,uiStopFrame,
            g_FGSegthresh,g_SegFoodThesMax,g_SegFoodThesMin,gthresEyeSeg,gthresEyeSegL,gi_MaxEllipseSamples,gi_VotesEllipseThres,gi_minEllipseMinor,gi_minEllipseMajor,gi_maxEllipseMajor,gi_CannyThresSmall,gi_CannyThres,
            gdMOGBGRatio,MOGhistory,gfVidfps
            ); // serialize things by passing them to the archive
     //archive();
     //archive(CEREAL_NVP(userROI));
   }

private:
   std::string strStateFilename = "trackerState.xml";
protected:
};

/// INIT FUNCTIONS ///


//void loadFromQrc(QString qrc,cv::Mat& imRes,int flag);

/// \brief Count Number of different Characters Between str1 and str2
int compString(QString str1,QString str2);

/// \brief sigsev Error handlers picking up unhandled errors
void installErrorHandlers();

/// \brief Process user provided config params and set corresponding internal/global variables
//void initGlobalParams(cv::CommandLineParser& parser,QStringList& inVidFileNames);

/// \brief Initializes ROI at start of tracking depending on user params / either large circle or user defined/configurable polygon
//void initROI();

/// \brief Initialize BG substractor objects, depending on options / can use cuda
void initBGSubstraction();

/// \brief Load internal and external template images memory cache //
int initDetectionTemplates();


float getAngleDiff(float anglefrom,float angleTo);
void testAngleDiff();

QString type2str(int type); //Opencv Mat type to string


#endif // CONFIG_H
