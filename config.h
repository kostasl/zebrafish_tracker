#ifndef CONFIG_H
#define CONFIG_H


#define ZTF_FISHCONTOURSIZE          80//40
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


#include <opencv2/core/core.hpp>
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


///

/// VIDEO AND BACKGROUND PROCESSING //
extern float gfVidfps;
extern const  unsigned int MOGhistory;//Use 100 frames Distributed across the video length To Find What the BGModel is
extern double gdMOGBGRatio; ///If a foreground pixel keeps semi-extern const ant value for about backgroundRatio*history frames, it's considered background and added to the model as a center of a new component.


//Processing Loop delay
extern uint cFrameDelayms;

extern const  double dLearningRate; //Learning Rate During Initial BG Modelling done over MOGhistory frames
extern const  double dLearningRateNominal;
extern double dBGMaskAccumulateSpeed;

/// BLOB DETECTION Filters //
//Area Filters
extern const  double dMeanBlobArea; //Initial Value that will get updated
extern const  double dVarBlobArea;
extern const  unsigned int gc_fishLength; //px Length Of Fish
extern const  unsigned int thresh_fishblobarea; //Min area above which to Filter The fish blobs
extern const  unsigned int thresh_maxfishblobarea; //Min area above which to Filter The fish blobs
extern const  unsigned int gthres_maxfoodblobarea;


/// extern const ants ///
extern int gcMaxFishModelInactiveFrames; //Number of frames inactive until track is deleted
extern int gcMaxFoodModelInactiveFrames; //Number of frames inactive (Not Matched to a Blob) until track is deleted
extern int gcMinFoodModelActiveFrames; //Min Number of consecutive frames it needs to be active  otherwise its deleted
extern int gi_FoodModelNumberLimit; //Maximum Number of Food Objects /Prey To track/Instantiate
extern const  int gMaxClusterRadiusFoodToBlob;
extern const  int thActive;// Deprecated If a track becomes inactive but it has been active less than thActive frames, the track will be deleted.
extern const  int thDistanceFish; //Threshold for distance between track-to blob assignement
//extern const  int thDistanceFood; //Threshold for distance between track-to blob assignement
extern int gFoodReportInterval;
extern const int gc_FishTailSpineSegmentLength_init;

extern const  int nTemplatesToLoad; //Number of Templates To Load Into Cache - These need to exist as images in QtResources


///Segmentation Params
extern int g_Segthresh; //Image Threshold to segment BG - Fish Segmentation uses a higher 2Xg_Segthresh threshold
extern int g_SegFoodThesMax; //Scan Thresholds to detect Food items in
extern int g_SegFoodThesMin;

extern int g_SegInnerthreshMult; //Image Threshold for Inner FIsh Features //Deprecated
extern int g_BGthresh; //BG threshold segmentation
extern int gi_ThresholdMatching; /// Minimum Score to accept that a contour has been found
extern bool gOptimizeShapeMatching; ///Set to false To disable matchShapes in FindMatching Contour
extern int gi_CannyThres;
extern int gi_CannyThresSmall; //Aperture size should be odd between 3 and 7 in function Canny
extern int gi_maxEllipseMajor; /// thres  for Eye Ellipse Detection methods
extern int gi_minEllipseMajor; ///thres for Eye Ellipse Detection methods (These Values Tested Worked Best)
extern int gi_minEllipseMinor; //The threshold in ellipse detection for how slim an ellipsoid can be
extern int gi_VotesEllipseThres; //Votes thres for The Backup Ellipse Detection Based on the Hough Transform
extern int gi_MaxEllipseSamples; //The number of fitted ellipsoids draw from the ranked queue to calculate mean fitted eye Ellipse
extern int gthresEyeSeg; //Central Threshold For Eye Segmentation In Isolated Head IMage
extern int gthresEyeSegL; //Lower Eye Seg Threshold
extern int gnumberOfTemplatesInCache; //INcreases As new Are Added
extern float gDisplacementThreshold; //Distance That Fish Is displaced so as to consider active and Record A point For the rendered Track /
extern int gFishBoundBoxSize; /// pixel width/radius of bounding Box When Isolating the fish's head From the image
extern int gFishTailSpineSegmentLength;
extern int gFitTailIntensityScanAngleDeg; //
extern double gUserReward; //User Given reward via GUI

extern const int gFishTailSpineSegmentCount;
extern const int gcFishContourSize;
extern const int gMaxFitIterations; //extern const ant For Max Iteration to Fit Tail Spine to Fish Contour

extern int giHeadIsolationMaskVOffset; //Vertical Distance to draw  Mask and Threshold Sampling Arc in Fish Head Mask
extern int giEyeIsolationMaskRadius; //Radius of circle masking between eyes
///Fish Features Detection Params
extern int gFishTemplateAngleSteps;
extern int gEyeTemplateAngleSteps;
extern double eyeStepIncrement;
extern double gTemplateMatchThreshold; //If not higher than 0.9 The fish body can be matched at extremeties
extern int iLastKnownGoodTemplateRow;
extern int iTemplateMatchFailCounter;
extern int iFishAngleOffset; //Manual User Corrected Fish Angle added to detected one
//using namespace std;

extern int iEyeMaskSepWidth; //Line width separating the eyes in the Head Image

///Contour Shaping - Global Variables
extern double sigma ;
extern int M ; //Gaussian Kernel Size

//Filename to Save Reinf. Learning State for EyeDetection
//extern QString gsEyeDetectorFilename;

/// \brief pointPairs is a vector holding points that the user chooses so as to measure distances between arbitrary points
typedef std::pair<cv::Point,cv::Point>  pointPair;
typedef std::vector<pointPair > pointPairs;


extern int keyboard; //input from keyboard
extern int screenx,screeny;
extern bool bshowMask; //True will show the BGSubstracted IMage/Processed Mask
extern bool bROIChanged;
extern bool bPaused;
extern bool bStartPaused;
extern bool bExiting;
extern bool bTracking;
extern bool bTrackFood;
extern bool bAddPreyManually ;
extern bool bMeasure2pDistance ; /// A mode allowing 2point distance measurement
extern bool bTrackFish;
extern bool bRecordToFile;
extern bool bSaveImages;
extern bool b1stPointSet;
extern bool bMouseLButtonDown;
//bool bSaveBlobsToFile; //Check in fnct processBlobs - saves output CSV
extern bool bEyesDetected ; ///Flip True to save eye shape feature for future detection
extern bool bStoreThisTemplate ;
extern bool bDraggingTemplateCentre ;
extern bool bUseEllipseEdgeFittingMethod  ; //Allow to Use the 2nd Efficient Method of Ellipsoid Fitting if the 1st one fails - Set to false to Make trakcing Faster
extern bool bFitSpineToTail ; // Runs The Contour And Tail Fitting Spine Optimization Algorith
extern bool bStartFrameChanged  ; /// When True, the Video Processing loop stops /and reloads video starting from new Start Position

extern bool bRenderToDisplay ; ///Updates Screen to User When True
extern bool bRenderWithAlpha; ///Blend Drawings so a A Transparency Effect is produced
extern bool bDrawFoodBlob    ; ///Draw circle around identified food blobs (prior to model matching)
extern bool bOffLineTracking ; ///Skip Frequent Display Updates So as to  Speed Up Tracking
extern bool bBlindSourceTracking ; /// Used for Data Labelling, so as to hide the data source/group/condition
extern bool bStaticAccumulatedBGMaskRemove ; /// Remove Pixs from FG mask that have been shown static in the Accumulated Mask after the BGLearning Phase
extern bool bUseBGModelling ; ///Use BG Modelling TO Segment FG Objects
extern bool gbUpdateBGModel ; //When Set a new BGModel Is learned at the beginning of the next video
extern bool gbUpdateBGModelOnAllVids; //When Set a new BGModel Is learned at the beginning of the next video
extern bool bApplyFishMaskBeforeFeatureDetection; ///Pass the masked image of the fish to the feature detector
extern bool bSkipExisting   ; /// If A Tracker DataFile Exists Then Skip This Video
extern bool bMakeCustomROIRegion; /// Uses Point array to construct
extern bool bUseMaskedFishForSpineDetect; /// When True, The Spine Is fit to the Masked Fish Image- Which Could Be problematic if The contour is not detected Well
extern bool bTemplateSearchThroughRows; /// Stops TemplateFind to Scan Through All Rows (diff temaplte images)- speeding up search + fail - Rows still Randomly Switch between attempts
extern bool bRemovePixelNoise; //Run Gaussian Filter Noise Reduction During Tracking
extern bool bUseGPU;
extern bool bUseOpenCL;
extern bool bUseHistEqualization; //To enhance to contrast in Eye Ellipse detection

extern int trackFnt;  //Font for Reporting - Tracking
extern float trackFntScale;

extern  std::string strTemplateImg; ///Load From Resource

extern  uint uiStartFrame;
extern  uint uiStopFrame;
//Contours Shaping
extern double sigma;
extern int M; //Gaussian Kernel Size


// Gaussian Curve Smoothing Kernels For fish Contour//
extern std::vector<double> gGaussian,dgGaussian,d2gGaussian;

extern QElapsedTimer gTimer;
extern QFile outfishdatafile;
extern QFile outfooddatafile;
extern QFile EyeDetectorRL; //Reinforcement Learned Behaviour For Eye Segmentation -

extern std::ofstream foutLog;//Used for Logging To File
static jmp_buf env; // For Skipping SIG sev and jumping to Prog Counter

extern QString outfilename;
extern std::string gstrwinName;
extern QString gstroutDirCSV,gstrinDirVid,gstrvidFilename; //The Output Directory


//Template size struct
extern cv::Size gszTemplateImg;
//For Inset Pasting
extern cv::Rect rect_pasteregion;

//Global Matrices Used to show debug images
extern cv::Mat frameDebugA,frameDebugB,frameDebugC,frameDebugD;
extern cv::Mat gframeCurrent,gframeLast; //Updated in processVideo Global Var Holding Copy of current and previous frame - usefull for opticflows

//Morphological Kernels
extern cv::Mat kernelOpen;
extern cv::Mat kernelDilateMOGMask;
extern cv::Mat kernelOpenfish;
extern cv::Mat kernelClose;
extern cv::Mat gLastfishimg_template;// OUr Fish Image Template
extern cv::Mat gFishTemplateCache; //A mosaic image contaning copies of template across different angles
//cv::Mat gEyeTemplateCache; //A mosaic image contaning copies of template across different angles


// Fish Detection CV BG Model//
extern cv::Ptr<cv::BackgroundSubtractorMOG2> pMOG2; //MOG2 Background subtractor
//cv::Ptr<cv::GeneralizedHough> pGHT;
//cv::Ptr<cv::GeneralizedHoughBallard> pGHTBallard;
//cv::Ptr<cv::GeneralizedHoughGuil> pGHTGuil;



/// \todo using a global var is a quick hack to transfer info from blob/Mask processing to fishmodel / Need to change the Blob Struct to do this properly
extern cv::Point gptHead; //Candidate Fish Contour Position Of HEad - Use for template Detect

extern ltROIlist vRoi;
//Rect Roi Keep Away from L-R Edges to Avoid Tracking IR lightRing Edges
extern  cv::Point ptROI1 ;
extern  cv::Point ptROI2;
extern cv::Point ptROI3 ;
extern cv::Point ptROI4 ;


class trackerState
{
 public:
 trackerState();
 void saveState(std::string strFilename);
 void loadState(std::string strFilename);


  int keyboard; //input from keyboard
  int screenx,screeny;
  // Option Flags //
  bool bshowMask; //True will show the BGSubstracted IMage/Processed Mask
  bool bStartPaused;
  bool bTracking;
  bool bTrackFood;
  bool bAddPreyManually ;
  bool bMeasure2pDistance ; /// A mode allowing 2point distance measurement
  bool bTrackFish;
  bool bRecordToFile;
  bool bSaveImages;
  bool bPaused;
  bool bUseEllipseEdgeFittingMethod  ; //Allow to Use the 2nd Efficient Method of Ellipsoid Fitting if the 1st one fails - Set to false to Make trakcing Faster
  bool bFitSpineToTail ; // Runs The Contour And Tail Fitting Spine Optimization Algorith
  bool bApplyFishMaskBeforeFeatureDetection; ///Pass the masked image of the fish to the feature detector
  bool bSkipExisting   ; /// If A Tracker DataFile Exists Then Skip This Video
  bool bMakeCustomROIRegion; /// Allows uset to set Custom ROI
  bool bUseMaskedFishForSpineDetect; /// When True, The Spine Is fit to the Masked Fish Image- Which Could Be problematic if The contour is not detected Well
  bool bTemplateSearchThroughRows; /// Stops TemplateFind to Scan Through All Rows (diff temaplte images)- speeding up search + fail - Rows still Randomly Switch between attempts
  bool bRemovePixelNoise; //Run Gaussian Filter Noise Reduction During Tracking
  bool bUseGPU;
  bool bUseOpenCL;
  bool bUseHistEqualization; //To enhance to contrast in Eye Ellipse detection
  std::string strTemplateImg; ///Load From Resource
  bool bBlindSourceTracking ; /// Used for Data Labelling, so as to hide the data source/group/condition
  bool bStaticAccumulatedBGMaskRemove ; /// Remove Pixs from FG mask that have been shown static in the Accumulated Mask after the BGLearning Phase
  bool bUseBGModelling ; ///Use BG Modelling TO Segment FG Objects
  bool gbUpdateBGModel ; //When Set a new BGModel Is learned at the beginning of the next video
  bool gbUpdateBGModelOnAllVids; //When Set a new BGModel Is learned at the beginning of the next video


  //Rect/Polygon Roi Keep Away from L-R Edges to Avoid Tracking IR lightRing Edges
  std::shared_ptr<std::vector<tROIpt>> userROI;

  std::string  gstroutDirCSV;//The Output Directory
  std::string  gstrinDirVid;
  std::string  gstrvidFilename;


  //Contours Shaping
  double curveSmoothKernelSigma;
  int curveSmoothKernelSize_M; //Gaussian Kernel Size


  //Screen Render Options//
  bool bRenderToDisplay ; ///Updates Screen to User When True
  bool bRenderWithAlpha;
  bool bOffLineTracking ; ///Skip Frequent Display Updates So as to  Speed Up Tracking
  bool bDrawFoodBlob    ; ///Draw circle around identified food blobs (prior to model matching)
  int trackFnt;  //Font for Reporting - Tracking
  float trackFntScale;


  ///Transient Action States //
  bool bMouseLButtonDown;
  bool b1stPointSet;
  bool bROIChanged;
  bool bExiting;
  bool bStoreThisTemplate ;
  bool bEyesDetected ; ///Flip True to save eye shape feature for future detection
  bool bDraggingTemplateCentre ;
  bool bStartFrameChanged  ; /// When True, the Video Processing loop stops /and reloads video starting from new Start Position
 //bool bSaveBlobsToFile; //Check in fnct processBlobs - saves output CSV


  ///Specific To the Tracked Video Options//
  uint uiStartFrame;
  uint uiStopFrame;
  QString outfilename;
  std::string gstrwinName;

    /// overload size operator / return full state object size
     size_t size() const _GLIBCXX_NOEXCEPT
    { return 1; //mStateValue.size()*mStateValue[1].size()*mStateValue[1][1].size();
     }

 // This method lets cereal know which data members to serialize
   template<class Archive>
   void serialize(Archive & archive)
   {
     archive(gstrwinName,CEREAL_NVP(gstroutDirCSV),CEREAL_NVP(gstrinDirVid),gstrvidFilename,
            bStartPaused,userROI,bRecordToFile,bTrackFish,bSaveImages,bUseEllipseEdgeFittingMethod,
            bTemplateSearchThroughRows,bApplyFishMaskBeforeFeatureDetection,bUseOpenCL,bUseGPU,bBlindSourceTracking,bStaticAccumulatedBGMaskRemove,
            gbUpdateBGModel,gbUpdateBGModelOnAllVids,bFitSpineToTail,bUseMaskedFishForSpineDetect,bUseHistEqualization,bRemovePixelNoise,bMeasure2pDistance,bAddPreyManually,
            bRenderToDisplay,bOffLineTracking,bDrawFoodBlob,curveSmoothKernelSigma,curveSmoothKernelSize_M
            ); // serialize things by passing them to the archive
     //archive();
     //archive(CEREAL_NVP(userROI));
   }

private:
   std::string strStateFilename = "trackerState.xml";
protected:
};

/// INIT FUNCTIONS ///

/// \brief Load Q Resources
void loadFromQrc(QString qrc,cv::Mat& imRes,int flag);

/// \brief Count Number of different Characters Between str1 and str2
int compString(QString str1,QString str2);

/// \brief sigsev Error handlers picking up unhandled errors
void installErrorHandlers();

/// \brief Process user provided config params and set corresponding internal/global variables
void initGlobalParams(cv::CommandLineParser& parser,QStringList& inVidFileNames);

/// \brief Initializes ROI at start of tracking depending on user params / either large circle or user defined/configurable polygon
void initROI();

/// \brief Initialize BG substractor objects, depending on options / can use cuda
void initBGSubstraction();

/// \brief Load internal and external template images memory cache //
int initDetectionTemplates();



#endif // CONFIG_H
