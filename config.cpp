#include <config.h>



/// VIDEO AND BACKGROUND PROCESSING //
float gfVidfps                  = 430;
const unsigned int MOGhistory   = 92;//Use 100 frames Distributed across the video length To Find What the BGModel is
double gdMOGBGRatio       = 0.95; ///If a foreground pixel keeps semi-constant value for about backgroundRatio*history frames, it's considered background and added to the model as a center of a new component.

//Processing Loop delay
uint cFrameDelayms              = 1;

const double dLearningRate                = 1.0/(MOGhistory); //Learning Rate During Initial BG Modelling done over MOGhistory frames
const double dLearningRateNominal         = 0.00001;


/// BLOB DETECTION Filters //
//Area Filters
const double dMeanBlobArea                    = 100; //Initial Value that will get updated
const double dVarBlobArea                     = 20;
const unsigned int gc_fishLength        = 100; //px Length Of Fish
const unsigned int thresh_fishblobarea  = 350; //Min area above which to Filter The fish blobs
const unsigned int thresh_maxfishblobarea = 2250; //Min area above which to Filter The fish blobs
const unsigned int gthres_maxfoodblobarea = 150;


/// Constants ///
const int gcMaxFishModelInactiveFrames  = 300; //Number of frames inactive until track is deleted
const int gcMaxFoodModelInactiveFrames  = 150; //Number of frames inactive (Not Matched to a Blob) until track is deleted
const int gcMinFoodModelActiveFrames    = 100; //Min Number of consecutive frames it needs to be active  otherwise its deleted
const int gMaxClusterRadiusFoodToBlob   = 4;
const int thActive                      = 0;// Deprecated If a track becomes inactive but it has been active less than thActive frames, the track will be deleted.
const int thDistanceFish                = 150; //Threshold for distance between track-to blob assignement
const int thDistanceFood                = 4; //Threshold for distance between track-to blob assignement
const int gFoodReportInterval           = (int)gfVidfps;

const int nTemplatesToLoad  = 19; //Number of Templates To Load Into Cache - These need to exist as images in QtResources


///Segmentation Params
int g_Segthresh             = 31; //Image Threshold to segment BG - Fish Segmentation uses a higher 2Xg_Segthresh threshold
int g_SegInnerthreshMult    = 3; //Image Threshold for Inner FIsh Features //Deprecated
int g_BGthresh              = 10; //BG threshold segmentation
int gi_ThresholdMatching    = 10; /// Minimum Score to accept that a contour has been found
bool gOptimizeShapeMatching = false; ///Set to false To disable matchShapes in FindMatching Contour
int gi_CannyThres           = 150;
int gi_CannyThresSmall      = 50; //Aperture size should be odd between 3 and 7 in function Canny
int gi_maxEllipseMajor      = 21; /// thres  for Eye Ellipse Detection methods
int gi_minEllipseMajor      = 12; ///thres for Eye Ellipse Detection methods (These Values Tested Worked Best)
int gi_VotesEllipseThres        = 9; //Votes thres for The Backup Ellipse Detection Based on the Hough Transform
int gthresEyeSeg                = 135; //Threshold For Eye Segmentation In Isolated Head IMage
int gnumberOfTemplatesInCache   = 0; //INcreases As new Are Added
float gDisplacementThreshold    = 2.0; //Distance That Fish Is displaced so as to consider active and Record A point For the rendered Track /
int gFishBoundBoxSize           = 24; /// pixel width/radius of bounding Box When Isolating the fish's head From the image
int gFishTailSpineSegmentLength     = 10;

int gFitTailIntensityScanAngleDeg   = 60; //
const int gFishTailSpineSegmentCount= ZTF_TAILSPINECOUNT;
const int gcFishContourSize         = ZTF_FISHCONTOURSIZE;
const int gMaxFitIterations         = ZTF_TAILFITMAXITERATIONS; //Constant For Max Iteration to Fit Tail Spine to Fish Contour

int giHeadIsolationMaskVOffset      = 8; //Vertical Distance to draw  Mask and Threshold Sampling Arc in Fish Head Mask

///Fish Features Detection Params
int gFishTemplateAngleSteps     = 1;
int gEyeTemplateAngleSteps      = 5;
double gTemplateMatchThreshold  = 0.89; //If not higher than 0.9 The fish body can be matched at extremeties
int iLastKnownGoodTemplateRow   = 0;
int iLastKnownGoodTemplateCol   = 0;
//using namespace std;

///Global Variables
 double sigma = 3.0;
 int M = round((3.0*sigma+1.0) / 2.0) * 2 - 1; //Gaussian Kernel Size
