#ifndef CONFIG_H
#define CONFIG_H

#include <limits>
#include <string>
#include <cmath>

#define ZTF_FISHCONTOURSIZE          55
#define ZTF_TAILFITMAXITERATIONS     10 //For Spine To Contour Tail Fitting
#define ZTF_TAILSPINECOUNT          8
#if defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    #define USE_CUDA
#else
    #undef USE_CUDA
#endif

#undef  USE_CUDA_FOR_TEMPLATE_MATCHING //Can Be Very Slow When Using Current Algorithm for searching through Cache
#define TEMPLATE_COL_SEARCH_REGION  15 //Optimization How Far to Search Along Col in the Cache From the Starting POint



/// VIDEO AND BACKGROUND PROCESSING //
extern float gfVidfps;
extern const  unsigned int MOGhistory;//Use 100 frames Distributed across the video length To Find What the BGModel is
extern double gdMOGBGRatio; ///If a foreground pixel keeps semi-extern const ant value for about backgroundRatio*history frames, it's considered background and added to the model as a center of a new component.

//Processing Loop delay
extern uint cFrameDelayms;

extern const  double dLearningRate; //Learning Rate During Initial BG Modelling done over MOGhistory frames
extern const  double dLearningRateNominal;


/// BLOB DETECTION Filters //
//Area Filters
extern const  double dMeanBlobArea; //Initial Value that will get updated
extern const  double dVarBlobArea;
extern const  unsigned int gc_fishLength; //px Length Of Fish
extern const  unsigned int thresh_fishblobarea; //Min area above which to Filter The fish blobs
extern const  unsigned int thresh_maxfishblobarea; //Min area above which to Filter The fish blobs
extern const  unsigned int gthres_maxfoodblobarea;


/// extern const ants ///
extern const  int gcMaxFishModelInactiveFrames; //Number of frames inactive until track is deleted
extern const  int gcMaxFoodModelInactiveFrames; //Number of frames inactive (Not Matched to a Blob) until track is deleted
extern const  int gcMinFoodModelActiveFrames; //Min Number of consecutive frames it needs to be active  otherwise its deleted
extern const  int gMaxClusterRadiusFoodToBlob;
extern const  int thActive;// Deprecated If a track becomes inactive but it has been active less than thActive frames, the track will be deleted.
extern const  int thDistanceFish; //Threshold for distance between track-to blob assignement
//extern const  int thDistanceFood; //Threshold for distance between track-to blob assignement
extern const  int gFoodReportInterval;

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
extern int gi_VotesEllipseThres; //Votes thres for The Backup Ellipse Detection Based on the Hough Transform
extern int gthresEyeSeg; //Threshold For Eye Segmentation In Isolated Head IMage
extern int gnumberOfTemplatesInCache; //INcreases As new Are Added
extern float gDisplacementThreshold; //Distance That Fish Is displaced so as to consider active and Record A point For the rendered Track /
extern int gFishBoundBoxSize; /// pixel width/radius of bounding Box When Isolating the fish's head From the image
extern int gFishTailSpineSegmentLength;
extern int gFitTailIntensityScanAngleDeg; //

extern const int gFishTailSpineSegmentCount;
extern const int gcFishContourSize;
extern const int gMaxFitIterations; //extern const ant For Max Iteration to Fit Tail Spine to Fish Contour

extern int giHeadIsolationMaskVOffset; //Vertical Distance to draw  Mask and Threshold Sampling Arc in Fish Head Mask

///Fish Features Detection Params
extern int gFishTemplateAngleSteps;
extern int gEyeTemplateAngleSteps;
extern double gTemplateMatchThreshold; //If not higher than 0.9 The fish body can be matched at extremeties
extern int iLastKnownGoodTemplateRow;
extern int iFishAngleOffset; //Manual User Corrected Fish Angle added to detected one
//using namespace std;

///Global Variables
extern double sigma ;
extern int M ; //Gaussian Kernel Size


#endif // CONFIG_H
