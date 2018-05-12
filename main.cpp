///*
//// 25/11/2015 : kostasl Testing OpenCV bg substraction - to process larva in vial recording timelapses.
 //// App uses BG substyraction MOG2, with a slow learning rate.
 //// then Uses Open and Close / Dilation contraction techniques to get rid of noise and fill gaps
 //// Then uses cvBlob library to track and count larva on BG-processed images.
 //// The lib sources have been included to the source and slightly modified in update tracks to fix a bug.
 ////
 ///* User:
 ///    * Chooses input video file, then on the second dialogue choose the text file to export track info in CSV format.
 ///    * The green box defines the region over which the larvae are counted-tracked and recorded to file.
 ///    * Once the video begins to show, use to left mouse clicks to define a new region in the image over which you want to count the larvae.
 ///    * Press p to pause Image. once paused:
 ///    * s to save snapshots in CSV outdir pics subfolder.
 ///    * 2 Left Clicks to define the 2 points of region-of interest for tracking.
 ///    * m to show the masked image of the larva against BG.
 ///    * t Start Tracking
 ///    * p to Pause
 ///    * r to UnPause/Run
 ///    * D to delete currently used template from cache
 ///    * T to save current tracked region as new template
 ///    * q Exit Quit application
 ///*
 ///* \note  Changing ROI hits SEG. FAULTs in update tracks of the library. So I made setting of ROI only once.
 ///* The Area is locked after t is pressed to start tracking. Still it fails even if I do it through cropping the images.
 ///* So I reverted to not tracking - as the code does not work well - I am recording blobs For now
 ///*
 ///*  Dependencies : opencv3
 ///*
 /// Added: Detection of stopped Larva or loss of features from BG Substraction - via mask correction
 ///    *Filter blobs and maintain separate lists for each class (food/fish)
 ///    * track blobs of different class (food/fish) separatelly so tracks do not interfere
 ///    *Issues:
 ///        *Multiple models for same blob
 ///    *Template Matching to spot fish across angles
 ///    *Ellipsoid fitting on edge points
 ///
 ///  Added Record of Food Count at regular intervals on each video in case
 ///          no fish is being tracked ROI - This should show the evolution of prey Count in time
 ///
 /// /note Found source of MultiProcess SegFault in OpenCL, Added try block on MOG2, and then flag to switch off OpenCL
 ///
 ////////


#include <larvatrack.h>
#include <ellipse_detect.h>
#include <template_detect.h>
#include <zfttracks.h>
#include <fgmaskprocessing.h>

#include <errorhandlers.h> // My Custom Mem Fault Handling Functions and Debug

#include <random>

#include <string.h>


#include <QDirIterator>
#include <QDir>
#include <QDebug>
//#include <QThread>
#include <QTime>

//Open CV

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
//#include <opencv2/bgsegm.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/video/background_segm.hpp>


#include <GUI/mainwindow.h>
///Curve Smoothing and Matching
#include <CSS/CurveCSS.h>

// Tracker Constant Defines
#include <config.h>


/// VIDEO AND BACKGROUND PROCESSING //
float gfVidfps                  = 430;
const unsigned int MOGhistory   = 100;//Use 100 frames Distributed across the video length To Find What the BGModel is

//Processing Loop delay
uint cFrameDelayms              = 1;

double dLearningRate                = 1.0/(2*MOGhistory); //Learning Rate During Initial BG Modelling done over MOGhistory frames
double dLearningRateNominal         = 0.00001;


/// BLOB DETECTION Filters //
//Area Filters
double dMeanBlobArea                    = 100; //Initial Value that will get updated
double dVarBlobArea                     = 20;
const unsigned int gc_fishLength        = 100; //px Length Of Fish
const unsigned int thresh_fishblobarea  = 350; //Min area above which to Filter The fish blobs
const unsigned int thresh_maxfishblobarea = 2250; //Min area above which to Filter The fish blobs
const unsigned int gthres_maxfoodblobarea = 150;


/// Constants ///
const int gcMaxFishModelInactiveFrames  = 300; //Number of frames inactive until track is deleted
const int gcMaxFoodModelInactiveFrames  = 50; //Number of frames inactive until track is deleted
const int gcMinFoodModelActiveFrames    = 10; //Min Number of consecutive frames it needs to be active otherwise its deleted
const int gMaxClusterRadiusFoodToBlob   = 8;
const int thActive                      = 0;// Deprecated If a track becomes inactive but it has been active less than thActive frames, the track will be deleted.
const int thDistanceFish                = 150; //Threshold for distance between track-to blob assignement
const int thDistanceFood                = 4; //Threshold for distance between track-to blob assignement
const int gFoodReportInterval           = (int)gfVidfps;

const int nTemplatesToLoad  = 19; //Number of Templates To Load Into Cache - These need to exist as images in QtResources


///Segmentation Params
int g_Segthresh             = 30; //Image Threshold to segment BG - Fish Segmentation uses a higher 2Xg_Segthresh threshold
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
int gFishBoundBoxSize           = 22; /// pixel width/radius of bounding Box When Isolating the fish's head From the image
int gFishTailSpineSegmentLength     = 9;
const int gFishTailSpineSegmentCount= ZTF_TAILSPINECOUNT;
int gFitTailIntensityScanAngleDeg   = 60; //

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
const double sigma = 3.0;
const int M = round((3.0*sigma+1.0) / 2.0) * 2 - 1; //Gaussian Kernel Size

 // Gaussian Curve Smoothing Kernels For fish Contour//
 std::vector<double> gGaussian,dgGaussian,d2gGaussian;



QElapsedTimer gTimer;
QFile outdatafile;
QString outfilename;
std::string gstrwinName = "FishFrame";
QString gstroutDirCSV,gstrvidFilename; //The Output Directory

//Global Matrices Used to show debug images
cv::Mat frameDebugA,frameDebugB,frameDebugC,frameDebugD;
cv::Size gszTemplateImg;

//cv::Ptr<cv::BackgroundSubtractor> pMOG; //MOG Background subtractor
cv::Ptr<cv::BackgroundSubtractorMOG2> pMOG2; //MOG2 Background subtractor
//cv::Ptr<cv::BackgroundSubtractorKNN> pKNN; //MOG Background subtractor
//cv::Ptr<cv::bgsegm::BackgroundSubtractorGMG> pGMG; //GMG Background subtractor

// Fish Detection //
Ptr<GeneralizedHough> pGHT;
Ptr<GeneralizedHoughBallard> pGHTBallard;
Ptr<GeneralizedHoughGuil> pGHTGuil;

//Morphological Kernels
cv::Mat kernelOpen;
cv::Mat kernelOpenLaplace;
cv::Mat kernelOpenfish;
cv::Mat kernelClose;
cv::Mat gLastfishimg_template;// OUr Fish Image Template
cv::Mat gFishTemplateCache; //A mosaic image contaning copies of template across different angles
cv::Mat gEyeTemplateCache; //A mosaic image contaning copies of template across different angles

/// \todo using a global var is a quick hack to transfer info from blob/Mask processing to fishmodel / Need to change the Blob Struct to do this properly
cv::Point gptHead; //Candidate Fish Contour Position Of HEad - Use for template Detect


ltROIlist vRoi;
//

cv::Point ptROI1 = cv::Point(gFishBoundBoxSize+1,82);
cv::Point ptROI2 = cv::Point(618,82);
cv::Point ptROI3 = cv::Point(618,441);
cv::Point ptROI4 = cv::Point(gFishBoundBoxSize+1,441);



//Structures to hold blobs & Tracks
//cvb::CvTracks fishtracks;
//cvb::CvTracks foodtracks;
//cvb::CvTracks tracks; ///All tracks

//The fish ones are then revaluated using simple thresholding to obtain more accurate contours
fishModels vfishmodels; //Vector containing live fish models
foodModels vfoodmodels;

uint gi_MaxFishID;
uint gi_MaxFoodID; //Declared in Model Header Files

MainWindow* pwindow_main = 0;

// Other fonts:
//   CV_FONT_HERSHEY_SIMPLEX, CV_FONT_HERSHEY_PLAIN,
//   CV_FONT_HERSHEY_DUPLEX, CV_FONT_HERSHEY_COMPLEX,
//   CV_FONT_HERSHEY_TRIPLEX, CV_FONT_HERSHEY_COMPLEX_SMALL,
//   CV_FONT_HERSHEY_SCRIPT_SIMPLEX, CV_FONT_HERSHEY_SCRIPT_COMPLEX
int trackFnt = CV_FONT_HERSHEY_SIMPLEX;  //Font for Reporting - Tracking
float trackFntScale = 0.6;

// Global Control Vars ///

int keyboard; //input from keyboard
int screenx,screeny;
bool bshowMask; //True will show the BGSubstracted IMage/Processed Mask
bool bROIChanged;
bool bPaused;
bool bExiting;
bool bTracking;
bool bTrackFood = true;
bool bSaveImages = false;
bool b1stPointSet;
bool bMouseLButtonDown;
bool bSaveBlobsToFile; //Check in fnct processBlobs - saves output CSV
bool bEyesDetected = false; ///Flip True to save eye shape feature for future detection
bool bStoreThisTemplate = false;
bool bDraggingTemplateCentre = false;
bool bUseEllipseEdgeFittingMethod =false; //Allow to Use the 2nd Efficient Method of Ellipsoid Fitting if the 1st one fails - Set to false to Make trakcing Faster
bool bFitSpineToTail = true; // Runs The Contour And Tail Fitting Spine Optimization Algorith
bool bStartFrameChanged = false; /// When True, the Video Processing loop stops /and reloads video starting from new Start Position
bool bRenderToDisplay = true; ///Updates Screen to User When True
bool bOffLineTracking = false; ///Skip Frequent Display Updates So as to  Speed Up Tracking
bool bStaticAccumulatedBGMaskRemove       = true; /// Remove Pixs from FG mask that have been shown static in the Accumulated Mask after the BGLearning Phase
bool gbUseBGModelling                     = true; ///Use BG Modelling TO Segment FG Objects
bool gbUpdateBGModel                      = true; //When Set a new BGModel Is learned at the beginning of the next video
bool gbUpdateBGModelOnAllVids             = true; //When Set a new BGModel Is learned at the beginning of the next video
bool bApplyFishMaskBeforeFeatureDetection = true; ///Pass the masked image of the fish to the feature detector
bool bSkipExisting                        = true; /// If A Tracker DataFile Exists Then Skip This Video
bool bMakeCustomROIRegion                 = false; /// Uses Point array to construct
bool bUseMaskedFishForSpineDetect         = true; /// When True, The Spine Is fit to the Masked Fish Image- Which Could Be problematic if The contour is not detected Well
bool bTemplateSearchThroughRows           = false; /// Stops TemplateFind to Scan Through All Rows (diff temaplte images)- speeding up search + fail - Rows still Randomly Switch between attempts
/// \todo Make this path relative or embed resource
//string strTemplateImg = "/home/kostasl/workspace/cam_preycapture/src/zebraprey_track/img/fishbody_tmp.pgm";
string strTemplateImg = ":/img/fishbody_tmp"; ///Load From Resource


void loadFromQrc(QString qrc,cv::Mat& imRes,int flag = IMREAD_COLOR)
{
    //double tic = double(getTickCount());

    QFile file(qrc);

    if(file.open(QIODevice::ReadOnly))
    {
        qint64 sz = file.size();
        std::vector<uchar> buf(sz);
        file.read((char*)buf.data(), sz);
        imRes = imdecode(buf, flag);
    }else
        std::cerr << " Could not load template file " << qrc.toStdString();

    //double toc = (double(getTickCount()) - tic) * 1000.0 / getTickFrequency();
    //qDebug() << "OpenCV loading time: " << toc;

}

/// MAIN FUNCTION - ENTRY POINT ////

jmp_buf env; //For Memory Exception Handling

//Count Number of different Characters Between str1 and str2
int compString(QString str1,QString str2)
{
    int ret =0;
  for (int j=0;j<std::min(str1.length(),str2.length());j++)
        if (str1.mid(j) != str2.mid(j))
            ret++;

  return ret;
}


int main(int argc, char *argv[])
{
    bROIChanged = false;
    bPaused = false;
    bshowMask = false;
    bTracking = true; //Start By Tracking by default
    bExiting    = false;

    std::ofstream foutLog;//Used for Logging To File
    // Get the rdbuf of clog.
    // We will need it to reset the value before exiting.
    auto old_rdbufclog = std::clog.rdbuf();
    auto old_rdbufcerr = std::cerr.rdbuf();


       // install Error/Seg Fault handler
    if (signal(SIGSEGV, handler) == SIG_ERR)
    {
        std::cerr << "**Error Setting SIGSEV simple handler! ::" << strsignal(SIGSEGV) << std::endl;
    }

    if (setjmp (env) == 0) {
      if (signal(SIGABRT, on_sigabrt) == SIG_ERR)
         {
             std::cerr << "**Error Setting SIGABRT simple handler! ::" << strsignal(SIGABRT) << std::endl;
         }
    }

    if (setjmp (env) == 0) {
      signal(SIGBUS, &handler);
      if (signal(SIGBUS, on_sigabrt) == SIG_ERR)
         {
             std::cerr << "**Error Setting SIGBUS simple handler! ::" << strsignal(SIGABRT) << std::endl;
         }
    }


    ///Install Error Hanlder //
    struct sigaction sigact;

     sigact.sa_sigaction = crit_err_hdlr;
     sigact.sa_flags = SA_RESTART | SA_SIGINFO;

     if (sigaction(SIGSEGV, &sigact, (struct sigaction *)NULL) != 0)
     {
      fprintf(stderr, "error setting signal handler for %d (%s)\n",
        SIGSEGV, strsignal(SIGSEGV));

      exit(EXIT_FAILURE);
     }

    /// ERROR HANDLER SIGSEV /////////


    QApplication app(argc, argv);
    //QQmlApplicationEngine engine;


    ///Parse Command line Args

    /// Handle Command Line Parameters //
    const cv::String keys =
        "{help h usage ? |    | print this help  message   }"
        "{outputdir   o |    | Dir where To save sequence of images }"
        "{invideofile v |    | Behavioural Video file to analyse }"
        "{invideolist f |    | A text file listing full path to video files to process }"
        "{startframe s | 1  | Video Will start by Skipping to this frame    }"
        "{stopframe p | 0  | Video Will stop at this frame    }"
        "{duration d | 0  | Number of frames to Track for starting from start frame }"
        "{logtofile l |    | Filename to save clog stream to }"
        "{ModelBG b | 1  | Learn and Substract Stationary Objects from Foreground mask}"
        "{SkipTracked t | 0  | Skip Previously Tracked Videos}"
        "{PolygonROI r | 0  | Use pointArray for Custom ROI Region}"
        "{ModelBGOnAllVids a | 1  | Only Update BGModel At start of vid when needed}"
        ;

    cv::CommandLineParser parser(argc, argv, keys);

    stringstream ssMsg;
    ssMsg<<"Zebrafish Behavioural Video Tracker"<< std::endl;
    ssMsg<<"--------------------------"<<std::endl;
    ssMsg<<"Author : Konstantinos Lagogiannis 2017"<<std::endl;
    ssMsg<<"./zebraprey_track <outfolder> <inVideoFile> <startframe=1> <stopframe=0> <duration=inf>"<<std::endl;
    ssMsg<<"(note: folder is automatically generated when absent)"<<std::endl;
    ssMsg << "Example: \n  Use checkFilesProcessed.sh script to generate list of videos to processes then execute as : " << std::endl;
    ssMsg << "./zebrafish_track -f=VidFilesToProcessSplit1.txt -o=/media/kostasl/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData/Tracked30-11-17/" << std::endl;

    parser.about(ssMsg.str() );

    if (parser.has("help"))
    {
        parser.printMessage();

        return 0;
    }

    /// End Of Parse COmmandLine//


    MainWindow window_main;
    pwindow_main = &window_main;
    window_main.show();

    QString outfilename;

    if (parser.has("outputdir"))
    {
        string soutFolder   = parser.get<string>("outputdir");
        std::clog << "Cmd Line OutDir : " << soutFolder <<std::endl;
        gstroutDirCSV  = QString::fromStdString(soutFolder);

    }
    else
    {
      outfilename  = QFileDialog::getSaveFileName(0, "Save tracks to output","VX_pos.csv", "CSV files (*.csv);", 0, 0); // getting the filename (full path)
      gstroutDirCSV = outfilename.left(outfilename.lastIndexOf("/"));
    }


    std::cout << "Csv Output Dir is " << gstroutDirCSV.toStdString()  << "\n " <<std::endl;


    QStringList inVidFileNames;
    if (parser.has("invideofile"))
    {
        inVidFileNames.append( QString::fromStdString(parser.get<string>("invideofile")) );
    }

    if (parser.has("invideolist"))
    {
        qDebug() << "Load Video File List " <<  QString::fromStdString(parser.get<string>("invideolist"));
        QFile fvidfile( QString::fromStdString(parser.get<string>("invideolist")) );
        if (fvidfile.exists())
        {
            fvidfile.open(QFile::ReadOnly);
            //QTextStream textStream(&fvidfile);
            while (!fvidfile.atEnd())
            {
                QString line = fvidfile.readLine().trimmed();
                if (line.isNull())
                    break;
                else
                    inVidFileNames.append(line);
            }
        }else
        {
            qWarning() << fvidfile.fileName() << " does not exist!";
        }
    }

    /// Setup Output Log File //
    if ( parser.has("logtofile") )
    {
        qDebug() << "Set Log File To " <<  QString::fromStdString( parser.get<string>("logtofile") );

        QFileInfo oLogPath( QString::fromStdString(parser.get<string>("logtofile") ) );
        if (!oLogPath.absoluteDir().exists())
            QDir().mkpath(oLogPath.absoluteDir().absolutePath()); //Make Path To Logs
        foutLog.open(oLogPath.absoluteFilePath().toStdString());

         // Set the rdbuf of clog.
         std::clog.rdbuf(foutLog.rdbuf());
         std::cerr.rdbuf(foutLog.rdbuf());
    }


    //Check If We Are BG Modelling / BEst to switch off when Labelling Hunting Events
    if (parser.has("ModelBG"))
        gbUseBGModelling = (parser.get<int>("ModelBG") == 1)?true:false;

    if (parser.has("ModelBGOnAllVids"))
        gbUpdateBGModelOnAllVids = (parser.get<int>("ModelBGOnAllVids") == 1)?true:false;

    if (parser.has("SkipTracked"))
        bSkipExisting = (parser.get<int>("SkipTracked") == 1)?true:false;

    if (parser.has("PolygonROI"))
        bMakeCustomROIRegion = (parser.get<int>("PolygonROI") == 1)?true:false;





    //If No video Files have been loaded then Give GUI to User
    if (inVidFileNames.empty())
            inVidFileNames =QFileDialog::getOpenFileNames(0, "Select videos to Process",gstroutDirCSV.toStdString().c_str(), "Video file (*.mpg *.avi *.mp4 *.h264 *.mkv *.tiff *.png *.jpg *.pgm)", 0, 0);


    // get the applications dir pah and expose it to QML
    //engine.load(QUrl(QStringLiteral("qrc:///main.qml")));

    gTimer.start();
    //create GUI windows


    //cv::namedWindow(gstrwinName,CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    //cv::namedWindow(gstrwinName + " FishOnly",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    //cv::namedWindow("Debug C",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);

    //cv::namedWindow("Ellipse fit",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    /// Histogram Window With Manual Threshold Setting for Eye Seg ///
    //cv::namedWindow("HeadHist",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    //cv::setMouseCallback("HeadHist", CallBackHistFunc, NULL);
    ////



#ifdef    _ZTFDEBUG_
    cv::namedWindow("Debug D",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    cv::namedWindow("Debug A",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    cv::namedWindow("Debug B",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);


    frameDebugA = cv::Mat::zeros(640, 480, CV_8U);
    frameDebugB = cv::Mat::zeros(640, 480, CV_8U);
    frameDebugD = cv::Mat::zeros(640, 480, CV_8U);
#endif

    frameDebugC = cv::Mat::zeros(640, 480, CV_8U);
    //set the callback function for any mouse event
    //cv::setMouseCallback(gstrwinName, CallBackFunc, NULL);


    /// Init Polygon ROI ///
    ///Make A Rectangular Roi Default //
    if (bMakeCustomROIRegion)
    {
        std::vector<cv::Point> vPolygon;
        vPolygon.push_back(ptROI1); vPolygon.push_back(ptROI2); vPolygon.push_back(ptROI3); vPolygon.push_back(ptROI4);
        ltROI rectRoi(vPolygon);
        vRoi.push_back(rectRoi);
     }

    /// create Gaussian Smoothing kernels for Contour //
    getGaussianDerivs(sigma,M,gGaussian,dgGaussian,d2gGaussian);


    /// create Background Subtractor objects
    //(int history=500, double varThreshold=16, bool detectShadows=true
    //OPENCV 3
    pMOG2 =  cv::createBackgroundSubtractorMOG2(MOGhistory, 16,false);
    pMOG2->setHistory(MOGhistory);
    pMOG2->setNMixtures(20);
    pMOG2->setBackgroundRatio(0.10);
    //pMOG2->setShadowValue(255);

    //double dmog2TG = pMOG2->getVarThresholdGen();
    //pMOG2->setVarThresholdGen(1.0);
    //double dmog2VT = pMOG2->getVarThreshold();
    //pMOG2->setVarThreshold(3.0);

    ///////////////////////////////////////
    /// Setup Fish Body Template Cache //

    int idxTempl;

    for (idxTempl=0; idxTempl<nTemplatesToLoad;idxTempl++)
    {
        loadFromQrc(QString::fromStdString(strTemplateImg + to_string(idxTempl+1) + std::string(".pgm")),gLastfishimg_template,IMREAD_GRAYSCALE); //  loadImage(strTemplateImg);
        if (gLastfishimg_template.empty())
        {
            std::cerr << "Could not load template" << std::endl;
            exit(-1);
        }

        addTemplateToCache(gLastfishimg_template,gFishTemplateCache,idxTempl); //Increments Index
    }
    gszTemplateImg.width = gLastfishimg_template.size().width; //Save TO Global Size Variable
    gszTemplateImg.height = gLastfishimg_template.size().height; //Save TO Global Size Variable


    int ifileCount = loadTemplatesFromDirectory(gstroutDirCSV + QString("/templates/"));
    pwindow_main->nFrame = 1;
    pwindow_main->LogEvent(QString::number(ifileCount+nTemplatesToLoad) + QString("# Templates Loaded "));

    /// END OF FISH TEMPLATES ///


    ///* Create Morphological Kernel Elements used in processFrame *///
    kernelOpen      = cv::getStructuringElement(cv::MORPH_CROSS,cv::Size(1,1),cv::Point(-1,-1));
    kernelOpenLaplace = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(1,1),cv::Point(-1,-1));
    kernelOpenfish  = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(3,3),cv::Point(-1,-1)); //Note When Using Grad Morp / and Low res images this needs to be 3,3
    kernelClose     = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(3,3),cv::Point(-1,-1));


    //unsigned int hWnd = cvGetWindowHandle(sgstrwinName);
//    try{ //Check If cv is compiled with QT support //Remove otherwise
//        cv::setWindowTitle(strwinName, outfilename.toStdString());

//     //trackVideofiles(window_main);

//    }catch(int e)
//    {
//        std::cerr << "OpenCV not compiled with QT support! can display overlay" <<std::endl;
//    }

    //trackImageSequencefiles(window_main);
    try{

        //app.exec();
        unsigned int uiStartFrame = parser.get<uint>("startframe");
        unsigned int uiStopFrame = parser.get<uint>("stopframe");
        std::clog << gTimer.elapsed()/60000.0 << " >>> Start frame: " << uiStartFrame << " StopFrame: " << uiStopFrame << " <<<<<<<<<"  << std::endl;
        trackVideofiles(window_main,gstroutDirCSV,inVidFileNames,uiStartFrame,uiStopFrame);

    }catch (char *e)
    {
        //printf("Exception Caught: %s\n",e);
        qDebug() << "[Error] >>> Exception Caught while processing: " << outdatafile.fileName();
        std::cerr << "[Error] Memory Allocation Error :" << e;
        //std::cerr << "Memory Allocation Error! - Exiting";
        std::cerr << "[Error] Close And Delete Current output file: " << outdatafile.fileName().toStdString() ;
        closeDataFile(outdatafile);
        removeDataFile(outdatafile);
        app.quit();

        std::exit(EXIT_FAILURE);
        return EXIT_FAILURE;
    }


    //destroy GUI windows

    //cv::waitKey(0);                                          // Wait for a keystroke in the window
   //pMOG2->getBackgroundImage();
    //pMOG->~BackgroundSubtractor();
    //pMOG2->~BackgroundSubtractor();
    //pKNN->~BackgroundSubtractor();
    //pGMG->~BackgroundSubtractor();

    //Empty The Track and blob vectors
    //cvb::cvReleaseTracks(tracks);
    //cvb::cvReleaseBlobs(blobs);



    std::cout << "Total processing time : mins " << gTimer.elapsed()/60000.0 << std::endl;
    std::clog << "Total processing time : mins " << gTimer.elapsed()/60000.0 << std::endl;
///Clean Up //

    frameDebugA.release();
    frameDebugB.release();
    frameDebugC.release();
    frameDebugD.release();



    ///* Create Morphological Kernel Elements used in processFrame *///
    kernelClose.release();
    kernelOpenfish.release();
    kernelOpenLaplace.release();
    kernelOpen.release();
    gLastfishimg_template.release();
    gEyeTemplateCache.release();
    gFishTemplateCache.release();


    //app.quit();
    window_main.close();
    cv::destroyAllWindows();


    // Reset the rdbuf of clog.
     std::clog.rdbuf(old_rdbufclog);
     std::cerr.rdbuf(old_rdbufcerr);


    app.quit();
    //Catch Any Mem Alloc Error
    std::exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;

}



unsigned int trackVideofiles(MainWindow& window_main,QString outputFileName,QStringList invideonames,unsigned int istartFrame = 0,unsigned int istopFrame = 0)
{
    cv::Mat fgMask;
    cv::Mat bgMask;

    QString invideoname = "*.mp4";
    QString nextvideoname;
    //Show Video list to process
    //std::cout << "Video List To process:" <<std::endl;
    window_main.LogEvent("Video List To process:");
    for (int i = 0; i<invideonames.size(); ++i)
    {
       invideoname = invideonames.at(i);
       //std::cout << "*" <<  invideoname.toStdString() << std::endl;
       window_main.LogEvent(invideoname );
    }


    //Go through Each Image/Video - Hold Last Frame N , make it the start of the next vid.
    for (int i = 0; i<invideonames.size() && !bExiting; ++i)
    {

       //Empty Vector of Fish Models - and Reset ID Counter // gi_MaxFoodID = gi_MaxFishID = 1; - Done In Release
       ReleaseFishModels(vfishmodels);
       ReleaseFoodModels(vfoodmodels);

       invideoname = invideonames.at(i);
       nextvideoname = invideonames.at(std::min(invideonames.length(),i+1));
       gstrvidFilename = invideoname; //Global

       std::clog << gTimer.elapsed()/60000.0 << " Now Processing : "<< invideoname.toStdString() << " StartFrame: " << istartFrame << std::endl;
       //cv::displayOverlay(gstrwinName,"file:" + invideoname.toStdString(), 10000 );

       ///Open Output File Check If We Skip Processed Files
       if ( !openDataFile(outputFileName,invideoname,outdatafile) )
            if (bSkipExisting) //Failed Due to Skip Flag
                 continue; //Do Next File




       // Removed If MOG Is not being Used Currently - Remember to Enable usage in enhanceMask if needed//
       if (gbUseBGModelling && gbUpdateBGModel)
       {
            getBGModelFromVideo(bgMask, window_main,invideoname,outfilename,MOGhistory);
            cv::bitwise_not ( bgMask, bgMask ); //Invert Accumulated MAsk TO Make it an Fg Mask



            //Next Video File Most Likely belongs to the same Experiment / So Do not Recalc the BG Model
            if (compString(invideoname,nextvideoname) < 3 && !gbUpdateBGModelOnAllVids)
                gbUpdateBGModel = false; //Turn Off BG Updates

       }

       QFileInfo fiVidFile(invideoname);
       window_main.setWindowTitle("Tracking:" + fiVidFile.completeBaseName() );
       window_main.nFrame = 0;
       window_main.tickProgress(); //Update Slider



       //if (bStaticAccumulatedBGMaskRemove) //Hide The PopUp
        //cv::destroyWindow("Accumulated BG Model");

        //Can Return 0 If DataFile Already Exists and bSkipExisting is true
        uint ret = processVideo(bgMask,window_main,invideoname,outdatafile,istartFrame,istopFrame);

        if (ret == 0)
            window_main.LogEvent(" [Error] Could not open Video file for last video");

        if (ret == 1)
        {
            if (!bSkipExisting)
                std::cerr << gTimer.elapsed()/60000.0 << " Error Occurred Could not open data file for last video" << std::endl;
            else
                window_main.LogEvent(" Skipping  previously processed Video."); // std::cerr << gTimer.elapsed()/60000.0 << " Error Occurred Could not process last video" << std::endl;

            continue; //Do Next File
        }
        istartFrame = 1; //Reset So Next Vid Starts From The Beginnning
        istopFrame = 0; //Rest So No Stopping On Next Video
    }
    return istartFrame;
}



void processFrame(MainWindow& window_main,const cv::Mat& frame,cv::Mat& bgMask, unsigned int nFrame,cv::Mat& outframe,cv::Mat& frameHead)
{
    cv::Mat frame_gray,fgFishMask,fgFishImgMasked;
    cv::Mat fgFoodMask;

    std::vector<cv::KeyPoint> ptFoodblobs;
    std::vector<cv::KeyPoint> ptFishblobs;

    std::vector<std::vector<cv::Point> > fishbodycontours;
    std::vector<cv::Vec4i> fishbodyhierarchy;

    unsigned int nLarva         =  0;
    //unsigned int nFood          =  0;
    double dblRatioPxChanged    =  0.0;

    QString frameNumberString;
    frameNumberString = QString::number(nFrame);

    //For Morphological Filter
    ////cv::Size sz = cv::Size(3,3);
    //frame.copyTo(inputframe); //Keep Original Before Painting anything on it
    //update the background model
    //OPEN CV 2.4
    // dLearningRate is now Nominal value
    /// Apply Histogram Equalization
      frame.copyTo(outframe); //Make Replicate On which we draw output


    ///DRAW ROI
    if (bRenderToDisplay)
        drawAllROI(outframe);


    //lplframe = frameMasked; //Convert to legacy format

    //cvb::CvBlobs blobs;
    ///DO Tracking
    if (bTracking)
    {
       //Simple Solution was to Use Contours To measure LUarvae
        //cvtColo frame_grey
        //Draw THe fish Masks more accuratelly by threshold detection - Enhances full fish body detection
    //    enhanceFishMask(outframe, fgMask,fishbodycontours,fishbodyhierarchy);// Add fish Blobs
        cv::cvtColor( frame, frame_gray, cv::COLOR_BGR2GRAY);


        // Update BG Substraction Model
        cv::Mat fgMask;
        pMOG2->apply(frame_gray,fgMask,dLearningRateNominal);

        //Combine Masks and Remove Stationary Learned Pixels From Mask
        if (bStaticAccumulatedBGMaskRemove)
        {
           //cv::bitwise_not(fgMask,fgMask);
            //gbMask Is Inverted Already So It Is The Accumulated FGMASK, and fgMask is the MOG Mask
           cv::bitwise_and(bgMask,fgMask,fgMask); //Only On Non Stationary pixels - Ie Fish Dissapears At boundary
        }

        enhanceMask(frame_gray,fgMask,fgFishMask,fgFoodMask,fishbodycontours, fishbodyhierarchy);
        //frameMasked = cv::Mat::zeros(frame.rows, frame.cols,CV_8UC3);

        if (bApplyFishMaskBeforeFeatureDetection)
            frame_gray.copyTo(fgFishImgMasked,fgFishMask); //fgFishMask //Use Enhanced Mask
        else
            frame_gray.copyTo(fgFishImgMasked); //fgFishMask //Use Enhanced Mask

        //outframe.copyTo(fgFoodImgMasked,fgFoodMask); //Use Enhanced Mask
        //show the current frame and the fg masks
        //cv::imshow(gstrwinName + " FishOnly",frameMasked);



        cv::Mat maskedImg_gray;
        /// Convert image to gray and blur it
        cv::cvtColor( frame, maskedImg_gray, cv::COLOR_BGR2GRAY );

       // Filters Blobs between fish and food - save into global vectors
        //processBlobs(&lplframe,fgMask, blobs,tracks,gstroutDirCSV,frameNumberString,dMeanBlobArea);

        //Can Use Fish Masked fgFishImgMasked - But Templates Dont Include The masking
        processFishBlobs(fgFishImgMasked,fgFishMask, outframe , ptFishblobs);
        nLarva = ptFishblobs.size();


        ///Update Fish Models Against Image and Tracks - Obtain Bearing Angle Using Template
        //Can Use Fish Masked - But Templates Dont Include The masking
        //UpdateFishModels(fgFishImgMasked,vfishmodels,ptFishblobs,nFrame,outframe);
        UpdateFishModels(maskedImg_gray,vfishmodels,ptFishblobs,nFrame,outframe);
        //If A fish Is Detected Then Draw Its tracks
        fishModels::iterator ft = vfishmodels.begin();

        while (ft != vfishmodels.end() && bRenderToDisplay) //Render All Fish
        {
            fishModel* pfish = ft->second;
            assert(pfish);
            zftRenderTrack(pfish->zTrack, frame, outframe,CV_TRACK_RENDER_PATH, CV_FONT_HERSHEY_PLAIN,trackFntScale+0.2 );
            ++ft;
        }

        ///\todo Keep A Global List of all tracks?

        /// Isolate Head, Get Eye models, and Get and draw Spine model
        if (nLarva > 0)
            //An Image Of the Full Fish Is best In this Case
            //Do Not Use Masked Fish Image For Spine Fitting
            detectZfishFeatures(window_main, frame_gray,outframe,frameHead,fgFishImgMasked,fishbodycontours,fishbodyhierarchy); //Creates & Updates Fish Models

        ///////  Process Food Blobs ////
        // Process Food blobs

        if (bTrackFood)
        {

            processFoodBlobs(fgFoodMask,fgFoodMask, outframe , ptFoodblobs); //Use Just The Mask
            UpdateFoodModels(maskedImg_gray,vfoodmodels,ptFoodblobs,nFrame,outframe);

            //If A fish Is Detected Then Draw Its tracks
            foodModels::iterator ft = vfoodmodels.begin();
            while (ft != vfoodmodels.end() && bRenderToDisplay)
            {
                foodModel* pfood = ft->second;
                assert(pfood);
                zftRenderTrack(pfood->zTrack, frame, outframe,CV_TRACK_RENDER_ID | CV_TRACK_RENDER_BOUNDING_BOX, CV_FONT_HERSHEY_PLAIN,trackFntScale );
                ++ft;
            }


        }


    } //If Tracking

    //fishbodycontours.clear();
    //fishbodyhierarchy.clear();
    //Save to Disk

    ///
    /// \brief drawFrameText
    if (bRenderToDisplay)
        drawFrameText(window_main,nFrame,nLarva,vfoodmodels.size(),outframe);

    if (bshowMask && bTracking)
        cv::imshow("Isolated Fish",fgFishImgMasked);


    fgFishImgMasked.release();
    fgFishMask.release();
    fgFishImgMasked.release();

   // int RefCount = frame_gray.u ? (frame_gray.u->refcount) : 0; //Its 1 at this point as required
    //assert(RefCount == 1);
    //qDebug() << "frame_gray.u->refcount:" << RefCount;

    frame_gray.release();


} //End Of Process Frame

///
/// \brief drawFrameText  ///TEXT INFO Put Info TextOn Frame
/// \param inFrame
/// \param frameNumberString
/// \param nLarva
/// \param nFood
/// \param outFrame
///
void drawFrameText(MainWindow& window_main, uint nFrame,uint nLarva,uint nFood,cv::Mat& outframe)
{

    //Frame Number
    std::stringstream ss;

    QString frameNumberString;
    frameNumberString = QString::number(nFrame);
    char buff[200];
    static double vm, rss;

    cv::rectangle(outframe, cv::Point(10, 2), cv::Point(100,20),
               CV_RGB(10,10,10), -1);
    cv::putText(outframe, frameNumberString.toStdString(),  cv::Point(15, 15),
            trackFnt, trackFntScale ,  CV_RGB(250,250,0));

    //Count on Original Frame
    std::stringstream strCount;
    strCount << "Nf:" << (nLarva) << " Nr:" << nFood;
    cv::rectangle(outframe, cv::Point(10, 25), cv::Point(80,45),  CV_RGB(10,10,10), -1);
    cv::putText(outframe, strCount.str(), cv::Point(15, 38),
           trackFnt, trackFntScale ,  CV_RGB(250,250,0));


    //Report Time
    std::sprintf(buff,"t: %0.2f",gTimer.elapsed()/(1000.0*60.0) );
    //strLearningRate << "dL:" << (double)(dLearningRate);
    cv::rectangle(outframe, cv::Point(10, 50), cv::Point(50,70), cv::Scalar(10,10,10), -1);
    cv::putText(outframe, buff, cv::Point(15, 63),
            trackFnt, trackFntScale , CV_RGB(250,250,0));

//Memory Reporting Moved to Main Loop
//    //Time Rate - conv from ms to minutes
//    ///Report Status  + Memory Usage
//    if ((nFrame%300) == 0 || nFrame == 1)
//    {
//        //Report
//        ss.str("");
//        std::clog << "Frame:" << nFrame << " Processing time (mins): " << gTimer.elapsed()/60000.0 << std::endl;
//        //THats In KiB units /So 1Million is A Gigabyte
//        std::clog << "#Fish " << nLarva << " #Food Blobs:" << nFood << std::endl;

//        ss << "Frame:" << nFrame << " Processing time (mins): " << gTimer.elapsed()/60000.0;
//        window_main.LogEvent(QString::fromStdString(ss.str()));

//        //Show Memory Consumption
//        ss.str("");
//        ss.precision(4);
//        process_mem_usage(vm, rss);
//        std::clog << "Memory VM: " << vm/1024.0 << "MB; RSS: " << rss/1024.0 << "MB" << std::endl;
//        ss  << "Memory VM: " << vm/1024.0 << "MB; RSS: " << rss/1024.0 << "MB";
//        window_main.LogEvent(QString::fromStdString(ss.str()));

//    }//Report on Next Frame
//    if ((nFrame%301) == 0 || nFrame == 2)
//    {
//        process_mem_usage(vm, rss);
//        std::clog << "Delta Memory VM: " << vm/1024.0 << "MB; RSS: " << rss/1024.0 << "MB" << std::endl;
//        //Show Memory Consumption
//        ss.str("");
//        ss.precision(4);
//        ss  << "D Memory VM: " << vm/1024.0 << "MB; RSS: " << rss/1024.0 << "MB";
//        window_main.LogEvent(QString::fromStdString(ss.str()));
//    }
//    ///Show Memory On Image Frame
//    std::sprintf(buff,"Vm: %0.2fMB;Rss:%0.2fMB",vm/1024.0,rss/1024.0);
//    cv::rectangle(outframe, cv::Point(5, 490), cv::Point(80,510), cv::Scalar(10,10,10), -1);
//    cv::putText(outframe, buff, cv::Point(10, 505),
//            trackFnt,trackFntScale , CV_RGB(10,250,0));
} //DrawFrameText



//
// Process Larva video, removing BG, detecting moving larva- Setting the learning rate will change the time required
// to remove a pupa from the scene -
//
unsigned int processVideo(cv::Mat& bgMask, MainWindow& window_main, QString videoFilename, QFile& outdatafile, unsigned int startFrameCount,unsigned int stopFrame=0)
{

    //Speed that stationary objects are removed
    cv::Mat frame,outframe,outframeHead;
    unsigned int nFrame = 0;
    unsigned int nErrorFrames = 0;
    outframeHead = cv::Mat::zeros(gszTemplateImg.height,gszTemplateImg.width,CV_8UC1); //Initiatialize to Avoid SegFaults
    bPaused =false; //Start Paused



    QString frameNumberString;
    //?Replicate FG Mask to method specific
    //fgMask.copyTo(fgMaskMOG2);
    //fgMask.copyTo(fgMaskMOG);
    //fgMask.copyTo(fgMaskGMG);


    //Make Variation of FileNames for other Output

    //QString trkoutFileCSV = outFileCSV;



    //create the capture object
    cv::VideoCapture capture(videoFilename.toStdString());
    if(!capture.isOpened())
    {
        //error in opening the video input
        std::cerr << gTimer.elapsed()/60000.0 << " [Error] Unable to open video file: " << videoFilename.toStdString() << std::endl;
        return 0;
        //std::exit(EXIT_FAILURE);
    }


    gfVidfps  = capture.get(CAP_PROP_FPS);
    uint totFrames = capture.get(CV_CAP_PROP_FRAME_COUNT);
    window_main.setTotalFrames(totFrames);
    window_main.nFrame = nFrame;
    window_main.LogEvent(" **Begin Processing: " + videoFilename);
    std::cout << " **Begin Processing: " << videoFilename.toStdString() << std::endl; //Show Vid Name To StdOUt
    window_main.stroutDirCSV = gstroutDirCSV;
    window_main.vidFilename = videoFilename;
    QString strMsg(  " Vid Fps:" + QString::number(gfVidfps) + " Total frames:" + QString::number(totFrames) + " Start:" + QString::number(startFrameCount));
    window_main.LogEvent(strMsg);

    //qDebug() << strMsg;


    // Open OutputFile

   // if (!openDataFile(trkoutFileCSV,videoFilename,outdatafile))
   //     return 1;

    outfilename = outdatafile.fileName();

    capture.set(CV_CAP_PROP_POS_FRAMES,startFrameCount);
    nFrame = capture.get(CV_CAP_PROP_POS_FRAMES);
    frameNumberString = QString::number(nFrame);

    //read input data. ESC or 'q' for quitting
    while( !bExiting && (char)keyboard != 27 )
    {

        /// Flow Control Code  - For When Looking at Specific Frame Region ///
        // 1st Check If user changed Frame - and go to that frame
        if (bStartFrameChanged)
        {
            nFrame = window_main.nFrame;
            capture.set(CV_CAP_PROP_POS_FRAMES,window_main.nFrame);
            bPaused = true;
            bTracking = bTracking; //Do Not Change
            //bStartFrameChanged = false; //This is Reset Once The frame Is captured
            //Since we are jumping Frames - The fish Models Are invalidated / Delete
            ReleaseFishModels(vfishmodels);
            ReleaseFoodModels(vfoodmodels);
        }

        if (!bPaused)
        {
            nFrame = capture.get(CV_CAP_PROP_POS_FRAMES);
            window_main.nFrame = nFrame; //Update The Frame Value Stored in Tracker Window
            window_main.tickProgress();
        }


        if (nFrame == startFrameCount)
        {
            bTracking = true;
        }

         frameNumberString = QString::number(nFrame); //Update Display String Holding FrameNumber

    if (!bPaused || bStartFrameChanged)
    {
        bStartFrameChanged = false; //Reset

        try //Try To Read The Image of that video Frame
        {
            //read the current frame
            if(!capture.read(frame))
            {
                if (nFrame == startFrameCount)
                {
                    std::cerr << gTimer.elapsed()/60000.0 << " " <<  nFrame << "# [Error]  Unable to read first frame." << std::endl;
                    nFrame = 0; //Signals To caller that video could not be loaded.
                    //Delete the Track File //
                    std::cerr << gTimer.elapsed()/60000.0 << " [Error] Problem with Tracking - Delete Data File To Signal its Not tracked" << std::endl;
                    removeDataFile(outdatafile);

                    exit(EXIT_FAILURE);
                }
                else //Not Stuck On 1st Frame / Maybe Vid Is Over?>
                {
                   std::cerr << gTimer.elapsed()/60000.0 << " [Error] " << nFrame << "# *Unable to read next frame." << std::endl;
                   std::clog << gTimer.elapsed()/60000.0 << " Reached " << nFrame << "# frame of " << totFrames <<  " of Video. Moving to next video." <<std::endl;
                   //assert(outframe.cols > 1);

                   if (nFrame < totFrames-1)
                   {
                       std::cerr << gTimer.elapsed()/60000.0 << " [Error] " << nFrame << " [Error] Stopped Tracking before End of Video - Delete Data File To Signal its Not tracked" << std::endl;
                       removeDataFile(outdatafile); //Delete The Output File
                   }
                   else
                   {
                       std::clog << gTimer.elapsed()/60000.0 << " [info] processVideo loop done on frame: " << nFrame << std::endl;
                         ::saveImage(frameNumberString,gstroutDirCSV,videoFilename,outframe);
                   }
                   //continue;
                   break;
               }

            } //Can't Read Next Frame
        }catch(const std::exception &e)
        {
            std::cerr << gTimer.elapsed()/60000.0 << " [Error] reading frame " << nFrame << " skipping." << std::endl;

            if (nFrame < totFrames)
                capture.set(CV_CAP_PROP_POS_FRAMES,nFrame+1);

            nErrorFrames++;
            if (nErrorFrames > 20) //Avoid Getting Stuck Here
            {
                // Too Many Error / Fail On Tracking
                std::cerr << gTimer.elapsed()/60000.0 << " [Error]  Problem with Tracking Too Many Read Frame Errors - Stopping Here and Deleting Data File To Signal Failure" << std::endl;
                removeDataFile(outdatafile);

                break;
            }
            else
                continue;
        }

    } //If Not Paused //

    //Check If StopFrame Reached And Pause
    if (nFrame == stopFrame && stopFrame > 0 && !bPaused)
    {
         bPaused = true; //Stop Here
         std::cout << nFrame << " Stop Frame Reached - Video Paused" <<std::endl;
         pwindow_main->LogEvent(QString(">>Stop Frame Reached - Video Paused<<"));
    }

        nErrorFrames = 0;


        //Make Global Roi on 1st frame if it doesn't prexist
        if (vRoi.size() == 0)
        {
            ptROI2.x = (frame.cols-gFishBoundBoxSize)/2;
            ptROI2.y = gszTemplateImg.height/3;
        //Add Global Roi - Center - Radius
            ltROI newROI(cv::Point(frame.cols/2,frame.rows/2),ptROI2);
            addROI(newROI);

            //Check If FG Mask Has Been Created - And Make A new One
            if (bgMask.cols == 0)
            {
                bgMask = cv::Mat::zeros(frame.rows,frame.cols,CV_8UC1);
                // Add Roi To Mask Otherwise Make On Based oN ROI
                cv::circle(bgMask,newROI.centre,newROI.radius,CV_RGB(255,255,255),-1);
            }
        }

        //Pass Processed bgMask which Is then passed on to enhanceMask
        processFrame(window_main,frame,bgMask,nFrame,outframe,outframeHead);
        if (bRenderToDisplay)
        {
            window_main.showVideoFrame(outframe,nFrame); //Show On QT Window
            window_main.showInsetimg(outframeHead);
        }

        // Switch Off Render To Display If In Offline Tracking Mode
        //-Placed Here to Allow 1Frame Short Toggles Of Rendering
        bRenderToDisplay = !bOffLineTracking;
        //frame.copyTo(frameDebugD);
        //cv::imshow("Debug D",frameDebugD);

        /// Report Memory Usage Periodically - Every realtime Second//
        if ((nFrame % (uint)gfVidfps) == 0 || nFrame == 2)
        {
            double rss,vm;
            process_mem_usage(vm, rss);
            //std::clog << "Delta Memory VM: " << vm/1024.0 << "MB; RSS: " << rss/1024.0 << "MB" << std::endl;
            //Show Memory Consumption

            std::stringstream ss;
            ss.precision(4);
            ss  << " [Progress] " << nFrame <<"/" << totFrames << " D Memory VM: " << vm/1024.0 << "MB; RSS: " << rss/1024.0 << "MB";
            window_main.LogEvent(QString::fromStdString(ss.str()));

            // Render Next Frame To Display
            if (bOffLineTracking)
                bRenderToDisplay = true;

        }

        if (bSaveImages)
        {
            cv::putText(frameDebugD, "REC", cv::Point(15, 20),
                    cv::FONT_HERSHEY_SIMPLEX, 0.8 , cv::Scalar(250,250,0));
            ::saveImage(frameNumberString,gstroutDirCSV,videoFilename,outframe);
        }


        if (bshowMask)
        {
           // cv::imshow(gstrwinName + " FG Mask", fgMask);
            //cv::imshow(gstrwinName + " FG Fish Mask", fgMaskFish);
        }
        // Show Debug Screens //
        //cv::imshow("Debug A",frameDebugA);
        //cv::imshow("Debug B",frameDebugB);
        //cv::imshow("Debug C",frameDebugC);



        //Save only when tracking - And Not While Paused
        if (bTracking && !bPaused)
            saveTracks(vfishmodels,outdatafile,frameNumberString);
            //saveTracks(vfishmodels,trkoutFileCSV,videoFilename,frameNumberString);


        checkPauseRun(&window_main,keyboard,nFrame);


    } //main While loop
    //delete capture object
    capture.release();



    std::clog << gTimer.elapsed()/60000.0 << "[Progress] Exiting video processing loop <<<" <<std::endl;
    //Dont Forget to Reset startFrameCount = 1 So Next Video Starts from Beginning
    //stopFrame       = 0;//No Stopping on NExt Video
    //Close File
    closeDataFile(outdatafile);

    return nFrame; //Return Number of Last Frame Processed
}




//Operator for Priority Ordering
bool operator<(const fishModel& a, const fishModel& b)
{
  return a.templateScore > b.templateScore; //Max Heap
}


/// \brief Defines search area region and runs template matching
/// pt Cant be closer to image border than gFishBoundBoxSize - If it is it will fixed to this distance
/// Returns Angle of Matched Template, centre of Detected Template and Match Score
double doTemplateMatchAroundPoint(const cv::Mat& maskedImg_gray,cv::Point pt,int& detectedAngle,cv::Point& detectedPoint ,cv::Mat& frameOut )
{
    /// Fix Bounds For Search point such that search temaplte region is not smaller than template size
    pt.x = (pt.x <= gLastfishimg_template.cols)?(gLastfishimg_template.cols+2): pt.x;
    pt.x = (maskedImg_gray.cols-pt.x <= gLastfishimg_template.cols)?maskedImg_gray.cols-gLastfishimg_template.cols-2: pt.x;

    pt.y = (pt.y <= gLastfishimg_template.rows)?(gLastfishimg_template.rows+2): pt.y;
    pt.y = (maskedImg_gray.rows-pt.y <= gLastfishimg_template.rows)?maskedImg_gray.rows-gLastfishimg_template.rows-2: pt.y;
    ///

    double maxMatchScore =0; //
    cv::Point gptmaxLoc; //point Of Bestr Match
    cv::Size szTempIcon(std::max(gLastfishimg_template.cols,gLastfishimg_template.rows),std::max(gLastfishimg_template.cols,gLastfishimg_template.rows));
    cv::Point rotCentre = cv::Point(szTempIcon.width/2,szTempIcon.height/2);
    ///
    /// Check If Track Centre Point Contains An image that matches a fish template
    /// \todo make HeadPoint/Tail point a Propery of FishBlob
    //cv::Point centroid = fishblob->pt;
     //Locate Centroid Region at a point between blob Centroid And Detect HeadPoint on Curve
   // cv::Point centroid = ((cv::Point)fishblob->pt-gptHead)/3+gptHead;
    cv::Point pBound1 = cv::Point(max(0,min(maskedImg_gray.cols,pt.x-gFishBoundBoxSize-2)), max(0,min(maskedImg_gray.rows,pt.y-gFishBoundBoxSize-2)));
    cv::Point pBound2 = cv::Point(max(0,min(maskedImg_gray.cols,pt.x+gFishBoundBoxSize+2)), max(0,min(maskedImg_gray.rows,pt.y+gFishBoundBoxSize+2)));

    // Look for Fish Template Within The Blob Region //
    cv::Rect rectFish(pBound1,pBound2);

    // Debug //
    //#ifdef _ZTFDEBUG_
    cv::rectangle(frameOut,rectFish,CV_RGB(20,200,150),1);
    //#endif

    cv::Mat fishRegion(maskedImg_gray,rectFish); //Get Sub Region Image


    //If blob exists but No Fish Model yet then Search Through Cache to improve matching;
    bool findBestMatch = (vfishmodels.size() == 0);
    if (findBestMatch)
        pwindow_main->LogEvent(QString("Look for Best Match in Templates"));

    int AngleIdx = templatefindFishInImage(fishRegion,gFishTemplateCache,szTempIcon, maxMatchScore, gptmaxLoc,iLastKnownGoodTemplateRow,iLastKnownGoodTemplateCol,findBestMatch);

    detectedAngle =AngleIdx*gFishTemplateAngleSteps;

    cv::Point top_left  = pBound1+gptmaxLoc; //Get top Left Corner Of Template Detected Region
    detectedPoint = top_left + rotCentre; //Get Centre Of Template Detection Region - Used for Tracking

    return maxMatchScore;
}


///
/// \brief UpdateFishModels starting from Blob Info do the processing steps to update FishModels for this frame,
/// \param maskedImg_gray
/// \param vfishmodels
/// \param fishblobs
/// \param nFrame
/// \param frameOut
///\todo - Add TimeOut Period Before Deleting Model
void UpdateFishModels(const cv::Mat& maskedImg_gray,fishModels& vfishmodels,zftblobs& fishblobs,unsigned int nFrame,cv::Mat& frameOut){

    qfishModels qfishrank;

    fishModel* pfish = NULL;

    fishModels::iterator ft;

     // Look through Blobs find Respective fish model attached or Create New Fish Model if missing
    for (zftblobs::iterator it = fishblobs.begin(); it!=fishblobs.end(); ++it)
    {
        zftblob* fishblob = &(*it);

        cv::Point ptbcentre = fishblob->pt; //Start As First Guess - This is updated When TemplMatching
        int bestAngle;
        double  maxMatchScore = 0.0;
        bool bModelFound = false;

        //Check Through Models And Find The Closest Fish To This FishBlob
        for ( ft  = vfishmodels.begin(); ft!=vfishmodels.end(); ++ft)
        {
             pfish = ft->second;
             ///Does this Blob Belong To A Known Fish Model?
             //Check Overlap Of This Model With The Blob - And Whether The Image of this Blob contains something That looks like a fish
             if (pfish->zfishBlob.overlap(pfish->zfishBlob,*fishblob) > 0 )
             {
                 //If Yes then assign the fish with the overlapping blob the template Match Score
                bModelFound = true;
                cv::Point ptSearch = pfish->ptRotCentre; //((cv::Point)fishblob->pt-gptHead)/3+gptHead;
                maxMatchScore = doTemplateMatchAroundPoint(maskedImg_gray,ptSearch,bestAngle,ptbcentre,frameOut);

                //Failed? Try the blob Head (From Enhance Mask) Detected position
                if ( maxMatchScore < gTemplateMatchThreshold)
                {
                  ptSearch = ((cv::Point)fishblob->pt-gptHead)/3+gptHead;
                  maxMatchScore = doTemplateMatchAroundPoint(maskedImg_gray,ptSearch,bestAngle,ptbcentre,frameOut);
                }
                pfish->templateScore = maxMatchScore;

                 if ( maxMatchScore >= gTemplateMatchThreshold)
                 {
                     //Some existing Fish Can be associated with this Blob - As it Overlaps from previous frame
                    ///Update Model State
                    // But not While it Is manually updating/ Modifying Bounding Box (Flags Are set in Mainwindow)
                    if (!bStoreThisTemplate && !bDraggingTemplateCentre) //Skip Updating Bound If this round we are saving The Updated Boundary
                    {
                        pfish->updateState(fishblob,maxMatchScore,bestAngle,ptbcentre,nFrame,gFishTailSpineSegmentLength,iLastKnownGoodTemplateRow,iLastKnownGoodTemplateCol);
                    }
                    else
                    { //Rotate Template Box - Since this cannot be done Manually
                        pfish->bearingAngle   = bestAngle;
                        pfish->bearingRads   =  bestAngle*CV_PI/180.0;
                    }

                 }
                 else //Below Thres Match Score
                 {
                       //Overide If We cant find that fish anymore/ Search from the start of the row across all angles
                       if (pfish->inactiveFrames > 3)
                           iLastKnownGoodTemplateCol = 0;
                         qDebug() << nFrame << " Guessing next TemplCol:" << iLastKnownGoodTemplateCol;
                 }

                 ////////  Write Angle / Show Box  //////
                 //Blobs may Overlap With Previously Found Fish But Match Score Is low - Then The Box Is still Drawn
                 pfish->drawBodyTemplateBounds(frameOut);
                //Add To Priority Q So we can Rank - Only If Blob Ovelaps ?
                qfishrank.push(pfish);
             }//if Models Blob Overlaps with this Blob


        } //For Each Fish Model

       //If the Blob Has no Model fish, and the template Match is low
       //then still create new model as this could be a fish we have not seen before -
       // And we avoid getting stuck searching for best model
       //
        if (!bModelFound) // && maxMatchScore >= gTemplateMatchThreshold  Model Does not exist for track - its a new track
        {
            //Make new fish Model
            //fishModel* fish= new fishModel(track,fishblob);
           fishModel* fish= new fishModel(*fishblob,bestAngle,ptbcentre);
           fish->ID = ++gi_MaxFishID;

           fish->updateState(fishblob,maxMatchScore,bestAngle,ptbcentre,nFrame,gFishTailSpineSegmentLength,iLastKnownGoodTemplateRow,iLastKnownGoodTemplateCol);

           vfishmodels.insert(IDFishModel(fish->ID,fish));
           qfishrank.push(fish); //Add To Priority Queue
           std::stringstream strmsg;
           strmsg << " New fishmodel: " << fish->ID << " with Template Score :" << fish->templateScore;
           //std::clog << nFrame << strmsg.str() << std::endl;
           pwindow_main->LogEvent(QString::fromStdString(strmsg.str()));

        }
//        //Report No Fish
        if (!bModelFound && maxMatchScore < gTemplateMatchThreshold )
        {
            std::clog << nFrame << "# Tscore:" << maxMatchScore << " No good match for Fish Found " << std::endl;

        }

    } //For Each Fish Blob

    ///\brief Check priority Queue Ranking Candidate Fish with TemplateSCore - Keep Top One Only
    fishModel* pfishBest = 0;
    double maxTemplateScore = 0.0;
    while (pfishBest==0 && qfishrank.size() > 0) //If Not In ROI Then Skip
    {
            pfishBest = qfishrank.top(); //Get Pointer To Best Scoring Fish
            ///Check If fish Model Is In ROI //
            for (std::vector<ltROI>::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
            {
                ltROI iroi = (ltROI)(*it);
                if (!iroi.contains(pfishBest->ptRotCentre))
                {
                    qfishrank.pop();
                    pfishBest =0;
                }
             }
   }//Search For Best Model

   if (pfishBest)
   {
        //qfishrank.pop();//Remove From Priority Queue Rank
        maxTemplateScore = pfishBest->templateScore;
        pfishBest->inactiveFrames   = 0; //Reset Counter
    }



    ///Delete All FishModels EXCEPT the best Match - Assume 1 Fish In scene / Always Retain 1 Model
    ft = vfishmodels.begin();
    while(ft != vfishmodels.end() ) //&& vfishmodels.size() > 1
    {
        pfish = ft->second;

        if (pfishBest != pfish && pfishBest != 0)
        {
            //Check Ranking Is OK, as long off course that a fishTemplate Has Been Found On This Round -
            //OtherWise Delete The model?
            //Assertion Fails When Old Model Goes Out Of scene and video Is retracked
            //assert(pfish->templateScore < maxTemplateScore || maxTemplateScore == 0);
            if (pfish->inactiveFrames > gcMaxFishModelInactiveFrames) //Check If it Timed Out / Then Delete
            {
                std::clog << gTimer.elapsed()/60000 << " " << nFrame << "# Deleted fishmodel: " << pfish->ID << " Low Template Score :" << pfish->templateScore << " when Best is :"<< maxTemplateScore << std::endl;
                ft = vfishmodels.erase(ft);
                delete(pfish);
                continue;
            }else
            {
                pfish->inactiveFrames ++; //Increment Time This Model Has Not Been Active
            }
        }
        ++ft; //Increment Iterator
    } //Loop To Delete Other FishModels



} //UpdateFishModels //




///
/// \brief UpdateFishModels starting from Blob Info do the processing steps to update FishModels for this frame,
/// \param maskedImg_gray
/// \param vfishmodels
/// \param fishblobs
/// \param nFrame
/// \param frameOut
///\todo - Add TimeOut Period Before Deleting Model
void UpdateFishModelsOrig(const cv::Mat& maskedImg_gray,fishModels& vfishmodels,zftblobs& fishblobs,unsigned int nFrame,cv::Mat& frameOut){

    qfishModels qfishrank;

    fishModel* pfish = NULL;

    fishModels::iterator ft;

    cv::Size szTempIcon(std::max(gLastfishimg_template.cols,gLastfishimg_template.rows),std::max(gLastfishimg_template.cols,gLastfishimg_template.rows));
    cv::Point rotCentre = cv::Point(szTempIcon.width/2,szTempIcon.height/2);

    cv::Point gptmaxLoc; //point Of Bestr Match

     // Look through Blobs find Respective fish model attached or Create New Fish Model if missing
    for (zftblobs::iterator it = fishblobs.begin(); it!=fishblobs.end(); ++it)
    {
        zftblob* fishblob = &(*it);
        ///
        /// Check If Track Centre Point Contains An image that matches a fish template
        /// \todo make HeadPoint/Tail point a Propery of FishBlob
        //cv::Point centroid = fishblob->pt;
         //Locate Centroid Region at a point between blob Centroid And Detect HeadPoint on Curve
        cv::Point centroid = ((cv::Point)fishblob->pt-gptHead)/3+gptHead;
        cv::Point pBound1 = cv::Point(max(0,min(maskedImg_gray.cols,centroid.x-gFishBoundBoxSize-2)), max(0,min(maskedImg_gray.rows,centroid.y-gFishBoundBoxSize-2)));
        cv::Point pBound2 = cv::Point(max(0,min(maskedImg_gray.cols,centroid.x+gFishBoundBoxSize+2)), max(0,min(maskedImg_gray.rows,centroid.y+gFishBoundBoxSize+2)));

        // Look for Fish Template Within The Blob Region //
        cv::Rect rectFish(pBound1,pBound2);

        // Debug //
//#ifdef _ZTFDEBUG_
        cv::rectangle(frameOut,rectFish,CV_RGB(20,200,150),1);
//#endif

        cv::Mat fishRegion(maskedImg_gray,rectFish); //Get Sub Region Image
        double maxMatchScore; //

        //If blob exists but No Fish Model yet then Search Through Cache to improve matching;
        bool findBestMatch = (vfishmodels.size() == 0);
        if (findBestMatch)
            pwindow_main->LogEvent(QString("Look for Best Match in Templates"));

        int AngleIdx = templatefindFishInImage(fishRegion,gFishTemplateCache,szTempIcon, maxMatchScore, gptmaxLoc,iLastKnownGoodTemplateRow,iLastKnownGoodTemplateCol,findBestMatch);

        int bestAngle =AngleIdx*gFishTemplateAngleSteps;
        cv::Point top_left = pBound1+gptmaxLoc;
        cv::Point ptbcentre = top_left + rotCentre;

        bool bModelFound = false;
        //Check Through Models And Find The Closest Fish To This FishBlob
        for ( ft  = vfishmodels.begin(); ft!=vfishmodels.end(); ++ft)
        {
             pfish = ft->second;
             ///Does this Blob Belong To A Known Fish Model?
             //Check Overlap Of This Model With The Blob - And Whether The Image of this Blob contains something That looks like a fish
             if (pfish->zfishBlob.overlap(pfish->zfishBlob,*fishblob) > 0 )
             {
                 //If Yes then assign the fish with the overlapping blob the template Match Score
                 bModelFound = true;
                 pfish->templateScore = maxMatchScore;
                 if ( maxMatchScore >= gTemplateMatchThreshold)
                 {
                     //Some existing Fish Can be associated with this Blob - As it Overlaps from previous frame
                    ///Update Model State
                    // But not While it Is manually updating/ Modifying Bounding Box (Flags Are set in Mainwindow)
                    if (!bStoreThisTemplate && !bDraggingTemplateCentre) //Skip Updating Bound If this round we are saving The Updated Boundary
                    {
                        pfish->updateState(fishblob,maxMatchScore,bestAngle,ptbcentre,nFrame,gFishTailSpineSegmentLength,iLastKnownGoodTemplateRow,iLastKnownGoodTemplateCol);
                    }
                    else
                    { //Rotate Template Box - Since this cannot be done Manually
                        pfish->bearingAngle   = bestAngle;
                        pfish->bearingRads   =  bestAngle*CV_PI/180.0;
                    }

                 }
                 else //Below Thres Match Score
                 {
                       //Overide If We cant find that fish anymore/ Search from the start of the row across all angles
                       if (pfish->inactiveFrames > 3)
                           iLastKnownGoodTemplateCol = 0;
                         qDebug() << nFrame << " Guessing next TemplCol:" << iLastKnownGoodTemplateCol;
                 }

                 ////////  Write Angle / Show Box  //////
                 //Blobs may Overlap With Previously Found Fish But Match Score Is low - Then The Box Is still Drawn
                 pfish->drawBodyTemplateBounds(frameOut);
                //Add To Priority Q So we can Rank - Only If Blob Ovelaps ?
                qfishrank.push(pfish);
             }//if Models Blob Overlaps with this Blob


        } //For Each Fish Model

       //If the Blob Has no Model fish, and the template Match is low
       //then still create new model as this could be a fish we have not seen before -
       // And we avoid getting stuck searching for best model
       //
        if (!bModelFound) // && maxMatchScore >= gTemplateMatchThreshold  Model Does not exist for track - its a new track
        {
            //Make new fish Model
            //fishModel* fish= new fishModel(track,fishblob);
           fishModel* fish= new fishModel(*fishblob,bestAngle,ptbcentre);
           fish->ID = ++gi_MaxFishID;

           fish->updateState(fishblob,maxMatchScore,bestAngle,ptbcentre,nFrame,gFishTailSpineSegmentLength,iLastKnownGoodTemplateRow,iLastKnownGoodTemplateCol);

           vfishmodels.insert(IDFishModel(fish->ID,fish));
           qfishrank.push(fish); //Add To Priority Queue
           std::stringstream strmsg;
           strmsg << " New fishmodel: " << fish->ID << " with Template Score :" << fish->templateScore;
           //std::clog << nFrame << strmsg.str() << std::endl;
           pwindow_main->LogEvent(QString::fromStdString(strmsg.str()));

        }
//        //Report No Fish
        if (!bModelFound && maxMatchScore < gTemplateMatchThreshold )
        {
            std::clog << nFrame << "# Tscore:" << maxMatchScore << " No good match for Fish Found " << std::endl;

        }

    } //For Each Fish Blob

    ///\brief Check priority Queue Ranking Candidate Fish with TemplateSCore - Keep Top One Only
    fishModel* pfishBest = 0;
    double maxTemplateScore = 0.0;
    while (pfishBest==0 && qfishrank.size() > 0) //If Not In ROI Then Skip
    {
            pfishBest = qfishrank.top(); //Get Pointer To Best Scoring Fish
            ///Check If fish Model Is In ROI //
            for (std::vector<ltROI>::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
            {
                ltROI iroi = (ltROI)(*it);
                if (!iroi.contains(pfishBest->ptRotCentre))
                {
                    qfishrank.pop();
                    pfishBest =0;
                }
             }
   }//Search For Best Model

   if (pfishBest)
   {
        //qfishrank.pop();//Remove From Priority Queue Rank
        maxTemplateScore = pfishBest->templateScore;
        pfishBest->inactiveFrames   = 0; //Reset Counter
    }



    ///Delete All FishModels EXCEPT the best Match - Assume 1 Fish In scene / Always Retain 1 Model
    ft = vfishmodels.begin();
    while(ft != vfishmodels.end() ) //&& vfishmodels.size() > 1
    {
        pfish = ft->second;

        if (pfishBest != pfish && pfishBest != 0)
        {
            //Check Ranking Is OK, as long off course that a fishTemplate Has Been Found On This Round -
            //OtherWise Delete The model?
            //Assertion Fails When Old Model Goes Out Of scene and video Is retracked
            //assert(pfish->templateScore < maxTemplateScore || maxTemplateScore == 0);
            if (pfish->inactiveFrames > gcMaxFishModelInactiveFrames) //Check If it Timed Out / Then Delete
            {
                std::clog << gTimer.elapsed()/60000 << " " << nFrame << "# Deleted fishmodel: " << pfish->ID << " Low Template Score :" << pfish->templateScore << " when Best is :"<< maxTemplateScore << std::endl;
                ft = vfishmodels.erase(ft);
                delete(pfish);
                continue;
            }else
            {
                pfish->inactiveFrames ++; //Increment Time This Model Has Not Been Active
            }
        }
        ++ft; //Increment Iterator
    } //Loop To Delete Other FishModels



} //UpdateFishModels //



///
/// \brief UpdateFoodModels
/// \param maskedImg_gray
/// \param vfoodmodels
/// \param foodblobs
/// \param nFrame
/// \param frameOut
void UpdateFoodModels(const cv::Mat& maskedImg_gray,foodModels& vfoodmodels,zfdblobs& foodblobs,unsigned int nFrame,cv::Mat& frameOut)
{
    qfoodModels qfoodrank;
    foodModel* pfood = NULL;

    foodModels::iterator ft;


    /// Assign Blobs To Food Models //
     // Look through Blobs find Respective food model attached or Create New Food Model if missing
    for (zfdblobs::iterator it = foodblobs.begin(); it!=foodblobs.end(); ++it)
    {
        zfdblob* foodblob = &(*it);

        cv::Point centroid = foodblob->pt;
        cv::Point pBound1 = cv::Point(max(0,min(maskedImg_gray.cols,centroid.x-5)), max(0,min(maskedImg_gray.rows,centroid.y-5)));
        cv::Point pBound2 = cv::Point(max(0,min(maskedImg_gray.cols,centroid.x+5)), max(0,min(maskedImg_gray.rows,centroid.y+5)));

        cv::Rect rectFood(pBound1,pBound2);
        cv::rectangle(frameOut,rectFood,CV_RGB(10,150,150),1);
        // Debug //
#ifdef _ZTFDEBUG_
        cv::Mat fishRegion(maskedImg_gray,rectFish); //Get Sub Region Image
#endif


        //Check Through Models And Find The Closest Food To This FoodBlob
        for ( ft  = vfoodmodels.begin(); ft!=vfoodmodels.end(); ++ft)
        {
             pfood = ft->second;
             bool bMatch = false;
             ///Does this Blob Belong To A Known Food Model?

             //Skip This food Model if it Has Already Been Assigned on this
             // Frame Unless We Paused And Stuck on the same Frame
             if ((nFrame - pfood->nLastUpdateFrame)==0 && bPaused == false)
                continue;

             pfood->blobMatchScore = 0;//Reset So We Can Rank this Match
             //Add points as each condition is met
             //Is it the same
             pfood->blobMatchScore += pfood->zfoodblob.size - foodblob->size;

            //Penalize no Overlap
             float overlap = pfood->zfoodblob.overlap(pfood->zfoodblob,*foodblob);
             pfood->blobMatchScore -=10.0*overlap;
             if (overlap > 0.0)
                    bMatch = true;

             //Cluster Blobs to one model if within a fixed Radius  That are close
             int fbdist = norm(pfood->zTrack.centroid-foodblob->pt);
             pfood->blobMatchScore +=fbdist;
             if (fbdist < gMaxClusterRadiusFoodToBlob)
                 bMatch = true;
             //Rank Up if this food model has been around for a while, instead of newly created one
             if (pfood->activeFrames > gcMinFoodModelActiveFrames )
                 pfood->blobMatchScore += pfood->activeFrames/gcMinFoodModelActiveFrames;


             //Consider only Food Models Within Region Of Blob
             if  (bMatch)
             {
                 qfoodrank.push(pfood);

             }
               //If Yes then assign the food with the overlapping blob the Match Score
               //Some existing food Can be associated with this Blob - As it Overlaps from previous frame
         } // Loop Through Food Models


        ///\brief Check priority Queue Ranking Candidate food with SCore - Keep Top One Only
        foodModel* pfoodBest = 0;
        if (qfoodrank.size() > 0)
        {
            pfoodBest = qfoodrank.top(); //Get Pointer To Best Scoring Fish
            //qrank.pop();//Remove From Priority Queue Rank
            pfoodBest->inactiveFrames   = 0; //Reset Counter
            pfoodBest->activeFrames = pfoodBest->activeFrames+1; //Increase Count Of Consecutive Active Frames
            pfoodBest->updateState(foodblob,0,pfoodBest->zfoodblob.pt,nFrame,pfoodBest->blobMatchScore);

        }else  ///No Food Model Found - Create A new One //
        {
            pfoodBest = new foodModel(*foodblob,++gi_MaxFoodID);

            vfoodmodels.insert(IDFoodModel(pfoodBest->ID,pfoodBest));
            std::stringstream strmsg;
            strmsg << "# New foodmodel: " << pfoodBest->ID << " N:" << vfoodmodels.size();
            std::clog << nFrame << strmsg.str() << std::endl;

            pfoodBest->updateState(foodblob,0,pfoodBest->zfoodblob.pt,nFrame,500);
        }

        clearpq2(qfoodrank);

    } // Loop Through BLOBS


    ///Delete All Inactive Food Models
    ft = vfoodmodels.begin();
    while(ft != vfoodmodels.end() ) //&& vfishmodels.size() > 1
    {
        pfood = ft->second;
        // Delete If Inactive For Too Long
        //Delete If Not Active for Long Enough between inactive periods / Track Unstable
        if (pfood->inactiveFrames > gcMaxFoodModelInactiveFrames ||
            (pfood->activeFrames < gcMinFoodModelActiveFrames && pfood->inactiveFrames > gcMaxFoodModelInactiveFrames)
            ) //Check If it Timed Out / Then Delete
        {
            ft = vfoodmodels.erase(ft);
            std::clog << nFrame << "# Delete foodmodel: " << pfood->ID << " N:" << vfoodmodels.size() << std::endl;
            delete(pfood);

            continue;
        }
        else //INcrease Inactive Frame Count For this Food Model
        {//If this Model Has not Been Used Here
            if (pfood->nLastUpdateFrame-nFrame > 1)
            {
                pfood->activeFrames = 0; //Reset Count Of Consecutive Active Frames
                pfood->inactiveFrames ++; //Increment Time This Model Has Not Been Active
            }
        }

        ++ft; //Increment Iterator
    } //Loop To Delete Inactive FoodModels



} //UpdateFoodModels //




void keyCommandFlag(MainWindow* win, int keyboard,unsigned int nFrame)
{

    //implemend Pause
    if ((char)keyboard == 'p')
    {
        //frame.copyTo(frameCpy);
        bPaused = true;

        std::cout << "Paused" << endl;
    }

    if ((char)keyboard == 'q')
    {
        bExiting = true;
        std::cout << "Quit" << endl;
    }

    //Make Frame rate faster
    if ((char)keyboard == '+')
        cFrameDelayms--;
    //Slower
    if ((char)keyboard == '-')
        cFrameDelayms++;


    if ((char)keyboard == 't') //Toggle Tracking
    {
        if (!bTracking)
        {
            iLastKnownGoodTemplateRow = 0; //Reset Row
            iLastKnownGoodTemplateCol = 0;
            pwindow_main->LogEvent(QString("Tracking ON"));
        }else
            pwindow_main->LogEvent(QString("Tracking OFF"));

        bTracking = !bTracking;
    }

    if ((char)keyboard == 'f') //Toggle FOOD Tracking
    {
        bTrackFood=!bTrackFood;

        if (bTrackFood)
            pwindow_main->LogEvent(QString("Track food ON"));
        else
            pwindow_main->LogEvent(QString("Track food OFF"));
    }

    if ((char)keyboard == '[') //Rotate Template AntiClock Wise
    {
        iLastKnownGoodTemplateCol--;
        iLastKnownGoodTemplateCol = max(0,iLastKnownGoodTemplateCol);
        pwindow_main->LogEvent(QString("User Rotated Template:")+QString::number(iLastKnownGoodTemplateCol)  );
    }

    if ((char)keyboard == ']') //Rotate Template ClockWise
    {
        iLastKnownGoodTemplateCol++;
        iLastKnownGoodTemplateCol = min(360/gFishTemplateAngleSteps , iLastKnownGoodTemplateCol);
        pwindow_main->LogEvent(QString("User Rotated Template:")+QString::number(iLastKnownGoodTemplateCol)  );
    }


    if ((char)keyboard == 's')
    {
        std::cout << "Save Image" << endl;
        bSaveImages = !bSaveImages;
        if (bSaveImages)
            win->saveScreenShot();

    }

    if ((char)keyboard == 'r')
    {
        std::cout << "Run" << endl;
        bPaused = false;
        bStartFrameChanged = false;
        gTimer.start();
    }

    if ((char)keyboard == 'R')
    {
             std::cout << "Reset Spines for All Fish Models-" << endl;
             for (fishModels::iterator it=vfishmodels.begin(); it!=vfishmodels.end(); ++it)
             {
                 fishModel* fish = (*it).second;
                   //Let ReleaseTracks Handle This
                  fish->resetSpine();
             }
             //ReleaseFishModels(vfishmodels);
    }


    if ((char)keyboard == 'd')
    {
      bOffLineTracking = !bOffLineTracking; //Main Loop Will handle this
      if (bOffLineTracking)
        pwindow_main->LogEvent(QString(">> Offline Tracking Mode ON <<"));
      else
        pwindow_main->LogEvent(QString("<< Offline Tracking Mode OFF >>"));
    }


    if ((char)keyboard == 'q')
        bExiting = true; //Main Loop Will handle this
         //break;


//    //if ((char)keyboard == 'c')
//    if (nFrame > 1)
//    {
//      //  cv::imshow(gstrwinName, frame);
//       win->showCVimg(frame); //Show On QT Window
//    }



    //Toggle Show the masked - where blob id really happens
    if ((char)keyboard == 'm')
    {
             std::cout << "Show Mask" << std::endl;
             bshowMask = !bshowMask;
    }

    ///Flip Save Feature On - This Will last only for a single frame
    if ((char)keyboard == 'e')
    {
             std::cout << "Save Eye Feature on next frame" << std::endl;
             bEyesDetected = true;
    }


    if ((char)keyboard == 'T')
    {
             std::cout << "Store next Image as Template" << std::endl;
             bStoreThisTemplate = !bStoreThisTemplate;
    }

    if ((char)keyboard == 'D')
    {
        bStoreThisTemplate = false;
        std::stringstream ss;
        ss << "Delete Currently Used Template Image idx:" << iLastKnownGoodTemplateRow;
        pwindow_main->LogEvent(QString::fromStdString(ss.str()));
        deleteTemplateRow(gLastfishimg_template,gFishTemplateCache,iLastKnownGoodTemplateRow);
        iLastKnownGoodTemplateRow = 0;
    }

    if ((char)keyboard == 'z')
    {
        bStoreThisTemplate = false;
        std::stringstream ss;
        int iNewTemplateRow = (rand() % static_cast<int>(gnumberOfTemplatesInCache - 0 + 1));//Start From RANDOM rOW On Next Search
        ss << "Reset Used Template idx:" << iLastKnownGoodTemplateRow << " to " << iNewTemplateRow;


        pwindow_main->LogEvent(QString::fromStdString(ss.str()));
        iLastKnownGoodTemplateRow = iNewTemplateRow;
        iLastKnownGoodTemplateCol = 0;
    }

}


void checkPauseRun(MainWindow* win, int keyboard,unsigned int nFrame)
{

//    int ms = 1;
//    struct timespec ts = { ms / 1000, (ms % 1000) * 1000 * 1000 };
//    nanosleep(&ts, NULL);
    ///Memory Crash Here ///
    ///
//    try
//    {
        QCoreApplication::processEvents(QEventLoop::AllEvents);
//    }catch(...)
//    {
        //std::cerr << "Event Processing Exception!" << std::endl;
//        qWarning() << "Event Processing Exception!";
//        win->LogEvent(QString("Event Processing Exception!"));

  //  }
       // cv::waitKey(1);

        //while (bPaused && !bExiting)
       // {


            //Wait Until Key to unpause is pressed
            //

    if (bPaused)
    { //Spend more time processing GUI events when Paused
        //keyboard = cv::waitKey( 1 );
        QTime dieTime= QTime::currentTime().addMSecs(20);
            while (QTime::currentTime() < dieTime)
                QCoreApplication::processEvents(QEventLoop::AllEvents);

        //keyCommandFlag(win,keyboard,nFrame);
    }           //


  //              cv::waitKey(100);


        //}

}

bool saveImage(QString frameNumberString,QString dirToSave,QString filenameVid,cv::Mat& img)
{
    cv::Mat image_to_write;
    //cv::cvtColor(img,image_to_write, cv::COLOR_RGB2BGR); //BGR For imWrite
    img.copyTo(image_to_write);
    //Make ROI dependent File Name
    QFileInfo fiVid(filenameVid);
    QString fileVidCoreName = fiVid.completeBaseName();

    //Save Output BG Masks
    //QString imageToSave =   QString::fromStdString( std::string("output_MOG_") + frameNumberString + std::string(".png"));

    dirToSave.append("/pics/");
    QString imageToSave =  fileVidCoreName + "_" + frameNumberString + ".png";
    imageToSave.prepend(dirToSave);

    if (!QDir(dirToSave).exists())
    {
        std::clog << "Make directory " << dirToSave.toStdString() << std::endl;
        QDir().mkpath(dirToSave);
    }

    bool saved = cv::imwrite(imageToSave.toStdString(), image_to_write);
    if(!saved) {
        cv::putText(img,"Failed to Save " + imageToSave.toStdString(), cv::Point(25, 25), cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(250,250,250));
        cv::putText(img,"Failed to Save" + imageToSave.toStdString(), cv::Point(25, 25), cv::FONT_HERSHEY_SIMPLEX, 0.4 , cv::Scalar(0,0,0));
        std::cerr << "Unable to save " << imageToSave.toStdString() << std::endl;
        pwindow_main->LogEvent(QString("Failed to Save Image File - Retry"));
        if (!cv::imwrite(imageToSave.toStdString(), image_to_write))
        {
            pwindow_main->LogEvent(QString("2nd Failed attempt to Save Image File"));
            return false;
        }
        else
            std::cout << "Saved image " << imageToSave.toStdString() <<std::endl;

    }
    else
    {
     std::cout << "Saved image " << imageToSave.toStdString() <<std::endl;
    }

    //cv::imshow("Saved Frame", img);

    image_to_write.deallocate();
    return true;
}


/// Updated Blob Processing
/// \brief processFishBlobs Finds blobs that belong to fish
/// \param frame
/// \param maskimg
/// \param frameOut //Output Image With FishBlob Rendered
/// \param ptFishblobs opencv keypoints vector of the Fish
/// \return
///
int processFishBlobs(cv::Mat& frame,cv::Mat& maskimg,cv::Mat& frameOut,std::vector<cv::KeyPoint>& ptFishblobs)
{

    std::vector<cv::KeyPoint> keypoints;
    //std::vector<cv::KeyPoint> keypoints_in_ROI;
    cv::SimpleBlobDetector::Params params;

    params.filterByCircularity  = false;
    params.filterByColor        = false;
    params.filterByConvexity    = false;

    //params.maxThreshold = 16;
    //params.minThreshold = 8;
    //params.thresholdStep = 2;

    // Filter by Area.
    params.filterByArea = true;
    params.minArea = thresh_fishblobarea/2.0;
    params.maxArea = thresh_maxfishblobarea;

    /////An inertia ratio of 0 will yield elongated blobs (closer to lines)
    ///  and an inertia ratio of 1 will yield blobs where the area is more concentrated toward the center (closer to circles).
    /// WARNING Enabling filterByInertia Causes A Crash - (in Close.s-> thread)
    params.filterByInertia      = false;
    params.maxInertiaRatio      = 0.6;
    params.minInertiaRatio      = 0.1;


    //params.filterByInertia = true;

    // Set up the detector with default parameters.
    cv::Ptr<cv::SimpleBlobDetector> detector = cv::SimpleBlobDetector::create(params);

    // Critical To Provide the Mask Image and not the full frame //
    detector->detect( maskimg, keypoints); //frameMask

    //Mask Is Ignored so Custom Solution Required
    //for (cv::KeyPoint &kp : keypoints)
    ptFishblobs.clear();
    for(int i=0;i<keypoints.size();i++)
    {
        cv::KeyPoint kp = keypoints[i];

        ///Go Through Each ROI and Render Blobs -
        unsigned int RoiID = 0;
        for (std::vector<ltROI>::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
        {
            ltROI iroi = (ltROI)(*it);
            RoiID++;
            //Keypoint is in ROI so Add To Masked
            if (iroi.contains(kp.pt))
                     ptFishblobs.push_back(kp);

            //int maskVal=(int)gframeMask.at<uchar>(kp.pt);
            //if (maskVal > 0)
             //keypoints_in_mask.push_back(kp);
        }
    }


    // Draw detected blobs as red circles.
    // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures the size of the circle corresponds to the size of blob
    //frame.copyTo(frameOut,maskimg); //mask Source Image
    cv::drawKeypoints( frameOut, ptFishblobs, frameOut, cv::Scalar(250,20,20), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS ); //cv::DrawMatchesFlags::DRAW_RICH_KEYPOINTS


    detector->clear();

}


/// Updated Blob Processing
/// \brief processFoodBlobs Finds blobs that belong to rotifers
/// \param frame
/// \param maskimg
/// \param frameOut //Output Image With FishBlob Rendered
/// \param ptFoodblobs opencv keypoints vector of the Fish
/// \return
///
int processFoodBlobs(const cv::Mat& frame,const cv::Mat& maskimg,cv::Mat& frameOut,std::vector<cv::KeyPoint>& ptFoodblobs)
{

    std::vector<cv::KeyPoint> keypoints;


    //std::vector<cv::KeyPoint> keypoints_in_ROI;
    cv::SimpleBlobDetector::Params params;

    //a circle has a circularity of 1,
    //circularity of a square is 0.785, and so on.

    params.filterByCircularity  = false;
    params.minCircularity       = 0.8;
    params.maxCircularity       = 1.0;

    params.filterByColor        = false;
    params.filterByConvexity    = false;

    //params.maxThreshold = 16;
    //params.minThreshold = 8;
    //params.thresholdStep = 2;

    // Filter by Area.
    params.filterByArea = true;
    params.minArea = 2;
    params.maxArea = gthres_maxfoodblobarea;

    /////An inertia ratio of 0 will yield elongated blobs (closer to lines)
    ///  and an inertia ratio of 1 will yield blobs where the area is more concentrated toward the center (closer to circles).
    params.filterByInertia      = false;
    params.maxInertiaRatio      = 1.0;
    params.minInertiaRatio      = 0.1;


    //params.filterByInertia = true;

    // Set up the detector with default parameters.
    cv::Ptr<cv::SimpleBlobDetector> detector = cv::SimpleBlobDetector::create(params);

    //\todo - Memory Crash Here - double free corruption

    assert(maskimg.depth() == CV_8U);
    detector->detect( maskimg, keypoints); //frameMask


    //Mask Is Ignored so Custom Solution Required
    //for (cv::KeyPoint &kp : keypoints)
    ptFoodblobs.clear();
    for(int i=0;i<keypoints.size();i++)
    {
        cv::KeyPoint kp = keypoints[i];

        ///Go Through Each ROI and Render Blobs - Split Between Fish and Food
        unsigned int RoiID = 0;
        for (std::vector<ltROI>::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
        {
            ltROI iroi = (ltROI)(*it);
            RoiID++;
            //Keypoint is in ROI so Add To Masked
            if (iroi.contains(kp.pt))
                     ptFoodblobs.push_back(kp);

            //int maskVal=(int)gframeMask.at<uchar>(kp.pt);
            //if (maskVal > 0)
             //keypoints_in_mask.push_back(kp);
        }
    }


    // Draw detected blobs as red circles.
    // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures the size of the circle corresponds to the size of blob
    //frame.copyTo(frameOut,maskimg); //mask Source Image
    cv::drawKeypoints( frameOut, ptFoodblobs, frameOut, cv::Scalar(0,120,200), cv::DrawMatchesFlags::DRAW_RICH_KEYPOINTS );


    detector->clear();

    return ptFoodblobs.size();

}

int saveTrackedBlobs(cvb::CvBlobs& blobs,QString filename,std::string frameNumber,ltROI& roi)
{
    int cnt = 0;
    int Vcnt = 1;
    bool bNewFileFlag = true;

    //Loop Over ROI
    Vcnt++; //Vial Count
    cnt = 0;

    QFile data(filename);
    if (data.exists())
        bNewFileFlag = false;

    if(data.open(QFile::WriteOnly |QFile::Append))
    {
        QTextStream output(&data);
        if (bNewFileFlag)
             output << "frameN,SerialN,BlobLabel,Centroid_X,Centroid_Y,Area\n" ;

        //Loop Over Blobs
        for (cvb::CvBlobs::const_iterator it=blobs.begin(); it!=blobs.end(); ++it)
        {

            cvb::CvBlob* cvB = it->second;
            cv::Point pnt;
            pnt.x = cvB->centroid.x;
            pnt.y = cvB->centroid.y;

            cnt++;

            if (roi.contains(pnt))
                //Printing the position information
                output << frameNumber.c_str() << "," << cnt <<","<< cvB->label << "," << cvB->centroid.x <<","<< cvB->centroid.y  <<","<< cvB->area  <<"\n";
          }


       data.close();

      }


    return cnt;
}

//Saves the total Number of Counted Blobs and Tracks only
int saveTrackedBlobsTotals(cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString filename,std::string frameNumber,ltROI& roi)
{

    bool bNewFileFlag = true;
    int cnt = 0;
    int Larvacnt = 0;
    cnt++;
    //cv::Rect iroi = (cv::Rect)(*it);

    QFile data(filename);
    if (data.exists())
        bNewFileFlag = false;

    if(data.open(QFile::WriteOnly |QFile::Append))
    {

        int blobCount = 0;
        int trackCount = 0;

        //Count Blobs in ROI
        for (cvb::CvBlobs::const_iterator it = blobs.begin(); it!=blobs.end(); ++it)
        {
            cvb::CvBlob* blob = it->second;
            cv::Point pnt;
            pnt.x = blob->centroid.x;
            pnt.y = blob->centroid.y;

            if (roi.contains(pnt))
                blobCount++;
        }

        //Count Tracks in ROI
        for (cvb::CvTracks::const_iterator it = tracks.begin(); it!=tracks.end(); ++it)
        {
            cvb::CvTrack* track = it->second;
            cv::Point pnt;
            pnt.x = track->centroid.x;
            pnt.y = track->centroid.y;

            if (roi.contains(pnt))
                trackCount++;
        }


        QTextStream output(&data);
        if (bNewFileFlag)
             output << "frameN,blobN,TracksN \n";

        output << frameNumber.c_str() << "," << blobCount << "," << trackCount <<"\n";
        Larvacnt +=blobCount;
        data.close();
    }


    return Larvacnt;
}


//std::vector<cvb::CvBlob*> getBlobsinROI(cvb::CvBlobs& blobs)
//{
    //std::vector<cvb::CvBlob*> *vfiltBlobs = new std::vector<cvb::CvBlob*>((blobs.size()));

   // return 0;

//}


ltROI* ltGetFirstROIContainingPoint(ltROIlist& vRoi ,cv::Point pnt)
{
    for (ltROIlist::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
    {
        ltROI* iroi = &(*it);
        if (iroi->contains(pnt))
                return iroi;
    }

    return 0; //Couldn't find it
}


bool openDataFile(QString filepathCSV,QString filenameVid,QFile& data)
{
    int Vcnt = 1;
    bool newFile = false;
    //Make ROI dependent File Name
    QFileInfo fiVid(filenameVid);
    QFileInfo fiOut(filepathCSV);
    QString fileVidCoreName = fiVid.completeBaseName();
    QString dirOutPath = fiOut.absolutePath() + "/"; //filenameCSV.left(filenameCSV.lastIndexOf("/")); //Get Output Directory

    char buff[50];
    sprintf(buff,"_tracks_%d.csv",Vcnt);
    dirOutPath.append(fileVidCoreName); //Append Vid Filename To Directory
    dirOutPath.append(buff); //Append extension track and ROI number

    data.setFileName(dirOutPath);
    //Make Sure We do not Overwrite existing Data Files
    while (!newFile)
    {
        if (!data.exists()) //Write HEader
        {
            newFile = true;

        }else{
            //File Exists
            if (bSkipExisting)
            {
                pwindow_main->LogEvent("[warning] Output File Exists and SkipExisting Mode is on.");
                std::cerr << "Skipping Previously Tracked Video File" << std::endl;
                return false; //File Exists Skip this Video

            }
            else
            {
                //- Create Name
            //FilenAme Is Linke AutoSet_12-10-17_WTNotFedRoti_154_002_tracks_1.csv
                //Increase Seq Number And Reconstruct Name
                Vcnt++;
                sprintf(buff,"_tracks_%d.csv",Vcnt);
                dirOutPath = fiOut.absolutePath() + "/";
                dirOutPath.append(fileVidCoreName); //Append Vid Filename To Directory
                dirOutPath.append(buff); //Append extension track and ROI number
                data.setFileName(dirOutPath);
            }
         }
    }
    if (!data.open(QFile::WriteOnly |QFile::Append))
    {

        std::cerr << "Could not open output file : " << std::endl;

        return false;
    }else {

        std::clog << "Opened file " << dirOutPath.toStdString() << " for data logging." << std::endl;

        /// Write Header //
        QTextStream output(&data);
        output << "frameN \t ROI \t fishID \t AngleDeg \t Centroid_X \t Centroid_Y \t EyeLDeg \t EyeRDeg \t ThetaSpine_0 \t ";
        for (int i=1;i<gFishTailSpineSegmentCount;i++)
            output <<  "DThetaSpine_" << i << "\t";

        output << " templateScore";
        output << "\t lastTailFitError";
        output << "\t lEyeFitScore";
        output << "\t rEyeFitScore";
        output << "\t nFailedEyeDetectionCount";
        output << "\t RotiferCount \n";
        //output.flush();

    }

    return true;
}

void closeDataFile(QFile& data)
{
    data.close();
    std::clog << gTimer.elapsed()/60000 << " Closed Output File " << data.fileName().toStdString() << std::endl;
}

void removeDataFile(QFile& data)
{

   std::clog << gTimer.elapsed()/60000 << "[Warning] Deleting Output File " << data.fileName().toStdString() << std::endl;

   if (data.exists())
    data.deleteLater();
}

///
/// \brief saveTracks -  record new fish track position - and rotifer count - only if fish is in view
/// \param vfish
/// \param data
/// \param frameNumber
/// \return
///
int saveTracks(fishModels& vfish,QFile& data,QString frameNumber)
{
    int cnt;
    int Vcnt = 0;

    //Loop Over ROI
    for (ltROIlist::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
    {
        cnt = 1;
        Vcnt++;
        ltROI iroi = (ltROI)(*it);
        //Make ROI dependent File Name

        QTextStream output(&data);

        //Save Tracks In ROI
        for (fishModels::iterator it=vfish.begin(); it!=vfish.end(); ++it)
        {
            cnt++;
            fishModel* pfish = it->second;
            cvb::CvLabel cvL = it->first;

            if (iroi.contains(pfish->ptRotCentre))
            {
                //Printing the position information +
                //+ lifetime; ///< Indicates how much frames the object has been in scene.
                //+active; ///< Indicates number of frames that has been active from last inactive period.
                //+ inactive; ///< Indicates number of frames that has been missing.
                output << frameNumber << "\t" << Vcnt  << "\t" << (*pfish);
                output << "\t" << vfoodmodels.size() << "\n";
            }
            //Empty Memory Once Logged
            pfish->zTrack.pointStack.clear();
            pfish->zTrack.pointStack.shrink_to_fit(); //Requires this Call of C++ otherwise It Doesnt clear
         }//For eAch Fish


         //If there is are no fish Then Add a regular Entry denoting the number of prey
        if (bTrackFood && vfish.size() == 0 && (frameNumber.toUInt()%gFoodReportInterval == 0 || frameNumber.toUInt()==1))
        {
            //make Null Fish
            fishModel* pNullfish = new fishModel();
            pNullfish->ID          = 0;
            pNullfish->bearingRads = 0.0;
            pNullfish->resetSpine();

            output << frameNumber << "\t" << Vcnt  << "\t" << (*pNullfish);
            output << "\t" << vfoodmodels.size() << "\n";
            delete pNullfish;
        }

   } //Loop ROI

 return cnt;
}


//Mouse Call Back Function
void CallBackFunc(int event, int x, int y, int flags, void* userdata)
{

    cv::Point ptMouse(x,y);

     if  ( event == cv::EVENT_LBUTTONDOWN )
     {
        bMouseLButtonDown = true;
         //ROI is locked once tracking begins
        ///CHANGE ROI Only when Paused and ONCE
        if (bPaused && !bROIChanged)
        { //Change 1st Point if not set or If 2nd one has been set
             if ( b1stPointSet == false)
             {
                ptROI1.x = x;
                ptROI1.y = y;
                //cv::circle(frameDebugA,ptROI1,3,cv::Scalar(255,0,0),1);

                b1stPointSet = true;
             }
             else //Second & Final Point
             {
                ptROI2.x = x;
                ptROI2.y = y;
                ltROI newROI(ptROI1,ptROI2);
                //roi = newROI;

                addROI(newROI);
                //drawROI(frame);
                b1stPointSet = false; //Rotate To 1st Point Again
             }
        }





        std::cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" <<std::endl;
     }

     if (event == cv::EVENT_LBUTTONUP)
     {
        bMouseLButtonDown = false;
     }
     else if  ( event == cv::EVENT_RBUTTONDOWN )
     {
        cv::Point mousepnt;
        mousepnt.x = x;
        mousepnt.y = y;
       std::cout << "Right button of the mouse is clicked - Delete ROI position (" << x << ", " << y << ")" <<std::endl;

        if (bPaused && !bROIChanged)
        {
            deleteROI(mousepnt);
            drawAllROI(frameDebugC);
        }
     }
     else if  ( event == cv::EVENT_MBUTTONDOWN )
     {
         std::cout << "Middle button of the mouse is clicked - position (" << x << ", " << y << ")" <<std::endl;
     }


     else if ( event == cv::EVENT_MOUSEMOVE )
     {

     }
}


//Mouse Call Back Function
void CallBackHistFunc(int event, int x, int y, int flags, void* userdata)
{

    if  ( event == cv::EVENT_LBUTTONUP )
     {
            cv::Point mousepnt;
            mousepnt.x = x;
            mousepnt.y = y;

            gthresEyeSeg = x;
            std::cout << "Eye Threshold Set to:" << gthresEyeSeg << std::endl;
    }
}

void addROI(ltROI& newRoi)
{
    //std::vector<cv::Rect>::iterator it= vRoi.end();
    //vRoi.insert(it,newRoi);
    vRoi.push_back(newRoi);
    //Draw the 2 points
    cv::circle(frameDebugC,ptROI1,3,cv::Scalar(255,0,0),1);
    cv::circle(frameDebugC,ptROI2,3,cv::Scalar(255,0,0),1);

   std::cout << "Added, total:" << vRoi.size() <<std::endl;

}

void deleteROI(cv::Point mousePos)
{
    std::vector<ltROI>::iterator it = vRoi.begin();

    while(it != vRoi.end())
    {
        ltROI* roi=&(*it);

        if (roi->contains(mousePos))
        {
            std::vector<ltROI>::iterator tmp = it;
            vRoi.erase(tmp);
           std::cout << "Deleted:" << roi->x() << " " << roi->y() <<std::endl;
            break;
        }
         ++it;

    }

}

void drawAllROI(cv::Mat& frame)
{
    //frameCpy.copyTo(frame); //Restore Original IMage
    for (std::vector<ltROI>::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
    {

        ltROI iroi = (ltROI)(*it);
         //cv::rectangle(frame,iroi,cv::Scalar(0,0,250));
        iroi.draw(frame);

        //Mark a centre to show that Tracking is ON / this ROI is being Tracked/Recorded
         if (bTracking)
         {
             cv::Point pt1;
             pt1.x = iroi.centre.x;
             pt1.y = iroi.centre.y;
             //pt2.x = pt1.x + iroi.radius;
             //pt2.y = pt1.y; //+ iroi.height;
             cv::circle(frame,pt1,3,cv::Scalar(255,0,0),1);
             //cv::circle(frame,pt2,3,cv::Scalar(255,0,0),1);

         }
    }
}


///
/// \brief findMatchingContour Looks for the inner contour in a 2 level hierarchy that matches the point coords
/// \param contours source array in which to search
/// \param hierarchy
/// \param pt - Position around which we are searching
/// \param level - The required hierarchy level description of the contour being searched for
/// \return Index of *child*/Leaf contour closest to point
///
int findMatchingContour(std::vector<std::vector<cv::Point> >& contours,
                              std::vector<cv::Vec4i>& hierarchy,
                              cv::Point pt,
                              int level)
{
    int idxContour           = -1;
    bool bContourfound       = false;
    int mindistToCentroid    = +10000; //Start Far
    int distToCentroid       = +10000;
    int matchContourDistance = 10000;

    ///Jump Over Strange Error
    if (hierarchy.size() !=contours.size())
        return -1;


    /// Render Only Countours that contain fish Blob centroid (Only Fish Countour)
   ///Search Through Contours - Draw contours + hull results

   ///Find Contour with Min Distance in shape and space -attach to closest contour
   //In Not found Search Again By distance tp Full Contour
       //Find Closest Contour
       for( int i = 0; i< (int)contours.size(); i++ )
       {

          //Filter According to desired Level
          if (level == 0) /////Only Process Parent Contours
          {
            if (hierarchy[i][3] != -1) // Need to have no parent
               continue;
            if (hierarchy[i][2] == -1)  // Need to have child
                continue;
            assert(hierarchy[hierarchy[i][2]][3] == i ); // check that the parent of the child is this contour i
          }

          if (level == 1) /////Only Process Child Contours
          {
              if (hierarchy[i][3] == -1) // Need to have a parent
                  continue;
//                   //Parent should be root
//                   if (hierarchy[hierarchy[i][3]][3] != -1)
//                       continue;
          }

          if (level == 2) ////Needs to be top Level Contour
          {
              if (hierarchy[i][3] != -1) // No Parent Contour
                  continue;
//                   //Parent should be root
//                   if (hierarchy[hierarchy[i][3]][3] != -1)
//                       continue;
          }



          //It returns positive (inside), negative (outside), or zero (on an edge)
          //Invert Sign and then Rank From Smallest to largest distance
          if (contours[i].size() > 0)
            matchContourDistance = distToCentroid = -cv::pointPolygonTest(contours[i],pt,true);

          //Measure Space Mod -Penalize Outside contour Hits - Convert Outside Contour Distances to X times further +ve (penalize)
          //Make Distance alway positive
          //matchContourDistance = (distToCentroid<0)?abs(distToCentroid)*20:distToCentroid;
          // qDebug() << "-c" << i << " D:" <<  distToCentroid;



          //Only Update if Spatial Distance is smaller but also moving from outside to inside of the shape
          //-ve is outside - 0 on border -
          //if(mindistToCentroid <= 0 && distToCentroid >= 0))
          {
               if (matchContourDistance < mindistToCentroid)
               {
                   //Otherwise Keep As blob Contour
                   idxContour = i;
                   mindistToCentroid = matchContourDistance;//New Min

                   //qDebug() << "-----MinTD:"<< matchContourDistance << "<- HDist:" << dHudist << " Sp:" << distToCentroid << "AreaDist:" << abs(tArea - sArea) << "LengthDist:" << abs(tLength - sLength);

                   //Reject match 0 in case contour is not actually there
                   //if (matchContourDistance < gi_ThresholdMatching)
                        bContourfound = true;
               }
           }
       }


   if (!bContourfound)
   {
       std::cerr << "Failed,Closest Contour :" << idxContour << " d:" << mindistToCentroid << std::endl;
       idxContour = -1;
   }
      //qDebug() << "-------Got best " <<  idxContour << " D:"<< mindistToCentroid;

   assert(idxContour < (int)contours.size());

   return idxContour;
}


///
/// \brief findMatchingContour Looks for the inner contour in a 2 level hierarchy that matches the point coords
/// \param contours source array in which to search
/// \param hierarchy
/// \param pt - Position around which we are searching
/// \param level - The required hierarchy level description of the contour being searched for
/// \param matchhull approx shape we are looking for
/// \param fittedEllipse Not Used - pointer to array of Rotated rect fitted ellipsoids
/// \return Index of *child*/Leaf contour closest to point
///
int findMatchingContour(std::vector<std::vector<cv::Point> >& contours,
                              std::vector<cv::Vec4i>& hierarchy,
                              cv::Point pt,
                              int level,
                              std::vector<cv::Point>& matchhull,
                              std::vector<cv::RotatedRect>& outfittedEllipse)
{
    int idxContour           = -1;
    bool bContourfound       = false;
    int mindistToCentroid    = +10000; //Start Far
    int distToCentroid       = +10000;
    int matchContourDistance = 10000;

    int tArea = 0;
    int sArea = 0;

    int tLength = 0;
    int sLength = 0;

    double dHudist = 0.0; //Shape Distance Hu moments distance measure from OpenCV

    /// Render Only Countours that contain fish Blob centroid (Only Fish Countour)
   ///Search Through Contours - Draw contours + hull results

   ///Find Contour with Min Distance in shape and space -attach to closest contour
   //In Not found Search Again By distance tp Full Contour
       //Find Closest Contour
       for( int i = 0; i< (int)contours.size(); i++ )
       {

          //Filter According to desired Level
          if (level == 0) /////Only Process Parent Contours
          {
            if (hierarchy[i][3] != -1) // Need to have no parent
               continue;
            if (hierarchy[i][2] == -1)  // Need to have child
                continue;
            assert(hierarchy[hierarchy[i][2]][3] == i ); // check that the parent of the child is this contour i
          }

          if (level == 1) /////Only Process Child Contours
          {
              if (hierarchy[i][3] == -1) // Need to have a parent
                  continue;
//                   //Parent should be root
//                   if (hierarchy[hierarchy[i][3]][3] != -1)
//                       continue;
          }

          if (level == 2) ////Needs to be top Level Contour
          {
              if (hierarchy[i][3] != -1) // No Parent Contour
                  continue;
//                   //Parent should be root
//                   if (hierarchy[hierarchy[i][3]][3] != -1)
//                       continue;
          }



          //It returns positive (inside), negative (outside), or zero (on an edge)
          //Invert Sign and then Rank From Smallest to largest distance
          if (contours[i].size() > 0)
            matchContourDistance = distToCentroid = -cv::pointPolygonTest(contours[i],pt,true);

          //Measure Space Mod -Penalize Outside contour Hits - Convert Outside Contour Distances to X times further +ve (penalize)
          //Make Distance alway positive
          //matchContourDistance = (distToCentroid<0)?abs(distToCentroid)*20:distToCentroid;
          // qDebug() << "-c" << i << " D:" <<  distToCentroid;


          ///Match Shape -
          /// \warning  If initial Shape Is not eye like this may be stuck into rejecting shapes
          // If A shape is provided
          //Find Contour Shape Similar to the one previously used for eye(ellipsoid)
          if (matchhull.size() > 5 && gOptimizeShapeMatching) //Only If Shape has been initialized/Given
          {
               dHudist = cv::matchShapes(matchhull,contours[i],CV_CONTOURS_MATCH_I2,0.0);
               matchContourDistance += dHudist*10.0; //Add Shape Distance onto / X Scale so it obtains relative importance
               // Now Check That distance is not too far otherwise reject shape
               //if (matchContourDistance > 1.0)
               //    continue; //Next Shape/Contour
               //qDebug() << "HuDist:" << dHudist*10.0;
             //Check Area

               tArea = cv::contourArea(matchhull);
               sArea = cv::contourArea(contours[i]);

               tLength = cv::arcLength(matchhull,true);
               sLength = cv::arcLength(contours[i],true);

               //Add Difference in Area to Score
               matchContourDistance += (abs(tArea - sArea));
              // qDebug() << "AreaDist:" << abs(tArea - sArea);

               matchContourDistance += abs(tLength - sLength);
               //qDebug() << "LengthDist:" << abs(tLength - sLength);

          }


          //Only Update if Spatial Distance is smaller but also moving from outside to inside of the shape
          //-ve is outside - 0 on border -
          //if(mindistToCentroid <= 0 && distToCentroid >= 0))
          {
               if (matchContourDistance < mindistToCentroid)
               {
                   //Otherwise Keep As blob Contour
                   idxContour = i;
                   mindistToCentroid = matchContourDistance;//New Min

                   //qDebug() << "-----MinTD:"<< matchContourDistance << "<- HDist:" << dHudist << " Sp:" << distToCentroid << "AreaDist:" << abs(tArea - sArea) << "LengthDist:" << abs(tLength - sLength);

                   //Reject match 0 in case contour is not actually there
                   //if (matchContourDistance < gi_ThresholdMatching)
                        bContourfound = true;
               }
           }
       }


   if (!bContourfound)
   {
       std::cerr << "Failed,Closest Contour :" << idxContour << " d:" << mindistToCentroid << std::endl;
       idxContour = -1;
   }
      //qDebug() << "-------Got best " <<  idxContour << " D:"<< mindistToCentroid;

   assert(idxContour < (int)contours.size());

   return idxContour;
}

///
/// \brief findIndexClosesttoPoint Returns Contour Index Closest To point pt
/// \param vPointChain
/// \param pt
/// \return Index of Vector Point closest to pt
///
int findIndexClosesttoPoint(std::vector<cv::Point> vPointChain,cv::Point pt)
{
    double dMindist = 10000.0;
    int iminIdx = 0;
    for (int i=0;i<vPointChain.size();i++)
    {
        double ddist = cv::norm(vPointChain[i]-pt);
        if (ddist < dMindist)
        {
            iminIdx = i;
            dMindist = ddist;
        }

    }

return iminIdx;
}

///
/// \brief enhanceFishMask Looks for fish countours and draws them onto the FG mask so as to enhance features
/// This is to recover Background substraction errors -
/// It then uses the fixed image to Find contours *main Internal and External fish contours* using on Masked Image Showing Fish Outline
/// \param frameImg - Raw Input camera input in Mat - colour or gray -
/// \param fgMask - Modified Enhanced FG Mask Image
/// \param outFishMask - Mask Enhanced for Fish Blob Detection
/// \param outFoodMask Enhanced for Food Blob Detection
/// \todo Cross Check Fish Contour With Model Position
/// - Tracker Picks Up Wrong contour Although Template Matching Finds the fish!
/// Note: Should Use MOG Mask for Blob Detect, But . But thresholded IMg For Countour FInding
void enhanceMask(const cv::Mat& frameImg, cv::Mat& fgMask,cv::Mat& outFishMask,cv::Mat& outFoodMask,std::vector<std::vector<cv::Point> >& outfishbodycontours, std::vector<cv::Vec4i>& outfishbodyhierarchy)
{

int max_thresh = 255;
cv::Mat frameImg_gray;
cv::Mat threshold_output;
cv::Mat threshold_output_COMB;
cv::Mat maskFGImg; //The FG Mask - After Combining Threshold Detection

//cv::imshow("MOG2 Mask Raw",maskFGImg);

///////////////// MOG Mask Is not Used Currently //
/////get rid of noise/food marks
////Apply Open Operation dilate(erode())
//cv::morphologyEx(maskFGImg,maskFGImg, cv::MORPH_OPEN, kernelOpen,cv::Point(-1,-1),1);
//////jOIN bLOB Do Close : erode(dilate())
//cv::morphologyEx(maskFGImg,maskFGImg, cv::MORPH_CLOSE, kernelClose,cv::Point(-1,-1),2);
/////Remove Speckles // Should leave fish INtact
//cv::filterSpeckles(maskFGImg,0,3,2 );
/////////// MOG Mask Is not Used Currently //


///// Convert image to gray, Mask and
//cv::cvtColor( frameImg, frameImg_gray, cv::COLOR_BGR2GRAY );
frameImg.copyTo(frameImg_gray); //Its Grey Anyway
//frameImg_gray = frameImg.clone();
//cv::GaussianBlur(frameImg_gray,frameImg_blur,cv::Size(3,3),0);

/// Detect edges using Threshold , A High And  low /
/// Trick, threshold Before Marking ROI - So as to Obtain Fish Features Outside Roi When Fish is incomplete Within The ROI
//cv::threshold( frameImg_gray, outFishMask, g_Segthresh*1.5, max_thresh, cv::THRESH_BINARY ); // Log Threshold Image + cv::THRESH_OTSU
outFishMask = cv::Mat::zeros(frameImg_gray.rows,frameImg_gray.cols,CV_8UC1);



cv::threshold( frameImg_gray, threshold_output, g_Segthresh, max_thresh, cv::THRESH_BINARY ); // Log Threshold Image + cv::THRESH_OTSU

//- Can Run Also Without THe BG Learning - But will detect imobile debri and noise MOG!
if (gbUseBGModelling && !fgMask.empty()) //We Have a (MOG) Model In fgMask - So Remove those Stationary Pixels
{
    cv::Mat fgMask_dilate; //Expand The MOG Mask And Intersect with Threshold
    //cv::morphologyEx(fgMask,fgMask_dilate,cv::MORPH_OPEN,kernelOpenfish,cv::Point(-1,-1),1);
    cv::dilate(fgMask,fgMask_dilate,kernelOpenfish,cv::Point(-1,-1),2);
    cv::bitwise_and(threshold_output,fgMask_dilate,maskFGImg); //Combine
    //fgMask.copyTo(maskFGImg);
}
else
{
    // Use thresh Only to Detect FG Fish Mask Detect
    threshold_output.copyTo(maskFGImg);
}

/// MASK FG ROI Region After Thresholding Masks - This Should Enforce ROI on Blob Detection  //
//frameImg_gray.copyTo(frameImg_gray,maskFGImg);


//cv::adaptiveThreshold(frameImg_gray, threshold_output,max_thresh,cv::ADAPTIVE_THRESH_MEAN_C,cv::THRESH_BINARY,g_Segthresh,0); //Last Param Is const substracted from mean
//ADAPTIVE_THRESH_MEAN_C

/////////////////Make Hollow Mask
//make Inner Fish MAsk /More Accurate Way
//cv::threshold( frameImg_gray, threshold_output_H, g_SegInnerthreshMult * g_Segthresh, max_thresh, cv::THRESH_BINARY ); //Log Threshold Image
//cv::dilate(threshold_output,threshold_output,kernelOpenfish,cv::Point(-1,-1),g_SegInnerthreshMult);
//Substract Inner from Outer
//cv::bitwise_xor(threshold_output,threshold_output_H,threshold_output_COMB);
///////////////////


//Make Hollow Mask Directly - Broad Approximate -> Grows outer boundary
//cv::dilate(maskFGImg,threshold_output,kernelOpenfish,cv::Point(-1,-1),1);
cv::morphologyEx(maskFGImg,threshold_output_COMB, cv::MORPH_GRADIENT, kernelOpenfish,cv::Point(-1,-1),1);

/// Find contours main Internal and External contour using on Masked Image Showing Fish Outline
/// //Used RETR_CCOMP that only considers 1 level children hierachy - I use the 1st child to obtain the body contour of the fish
outfishbodycontours.clear();
std::vector<std::vector<cv::Point> > fishbodycontours;
std::vector<cv::Vec4i> fishbodyhierarchy;

//Then Use ThresholdImage TO Trace More detailed Contours
//cv::dilate(threshold_output_COMB,threshold_output_COMB_fish,kernelOpenfish,cv::Point(-1,-1),4);
cv::findContours( threshold_output_COMB, fishbodycontours,fishbodyhierarchy, cv::RETR_CCOMP,cv::CHAIN_APPROX_SIMPLE , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE

//Make Food Mask OUt Of FG Model /After Removing Noise
//cv::erode(maskFGImg,maskFGImg,kernelClose,cv::Point(-1,-1),1);
cv::morphologyEx(maskFGImg,outFoodMask,cv::MORPH_OPEN,kernelOpen,cv::Point(-1,-1),1);
//cv::dilate(outFoodMask,outFoodMask,kernelClose,cv::Point(-1,-1),1);

//threshold_output_COMB.copyTo(outFoodMask);
//outFoodMask = maskFGImg.clone();

//cv::bitwise_and(outFoodMask,maskFGImg,outFoodMask); //Remove Regions OUtside ROI
//std::vector< std::vector<cv::Point> > fishbodyContour_smooth;

///Draw Only the largest contours that should belong to fish
/// \todo Other Match Shapes Could be used here
/// \todo Use WaterShed - Let MOG mask Be FG label and then watershed
//int idxFishContour = -1;
std::vector<cv::Point> curve; // THe Fish Contour to use for new Mask
for (int kk=0; kk< (int)fishbodycontours.size();kk++)
{
    curve.clear();

        ///Filter for what looks like a fish //
        /// Can use many methods here such as match shapes / Hashing etc.

        //Find Parent Contour
        if (fishbodyhierarchy[kk][3] != -1) // Need to have no parent
           continue;
        if (fishbodyhierarchy[kk][2] == -1)  // Need to have child
            continue;

        /// Lets try simple area filter - Assume no large object need to be BG substracted
        uint area  = cv::contourArea(fishbodycontours[kk]);



        cv::Point centroid;
        ///Check Area and then  Find the thresholded Fish Contour std::max(dMeanBlobArea*8,(double)thresh_fishblobarea)
        if (area >  thresh_fishblobarea && area <= thresh_maxfishblobarea) //If Contour Is large Enough then Must be fish
        {
            cv::Moments moments =  cv::moments(fishbodycontours[kk]);

            centroid.x = moments.m10/moments.m00;
            centroid.y = moments.m01/moments.m00;
            //If Contained In ROI
            for (std::vector<ltROI>::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
            {
                ltROI iroi = (ltROI)(*it);
                //Keypoint is in ROI so Add To Masked
                if (iroi.contains(centroid))
                {
                     curve = fishbodycontours[kk];
                     outfishbodyhierarchy.push_back(fishbodyhierarchy[kk]); //Save Hierarchy Too
                }
            }

//            curve = fishbodycontours[kk];
            //std::vector<cv::RotatedRect> rectFeatures;
            //Add Blob To candidate Region of interest Mask
            //idxFishContour = findMatchingContour(fishbodycontours,fishbodyhierarchy,centroid,-1,fgMaskcontours[kk],rectFeatures);
            //std::clog << " cntr area :" << area << std::endl;
        }

        //Skip Very Small Curves
        if (curve.size() < gcFishContourSize)
            continue;

        assert(M % 2 == 1); //M is an odd number

        ///// SMOOTH COntours /////

        vector<double> curvex,curvey,smoothx,smoothy,resampledcurveX,resampledcurveY ;
        PolyLineSplit(curve,curvex,curvey);

        std::vector<double> X,XX,Y,YY;
        std::vector<double> dXY;

        getdXcurve(curvex,sigma,smoothx,X,XX,gGaussian,dgGaussian,d2gGaussian,false);
        getdXcurve(curvey,sigma,smoothy,Y,YY,gGaussian,dgGaussian,d2gGaussian,false);

        dXY.resize(X.size());
        /// Find Tail As POint Of Maximum Curvature dXY
        int idxMin2 = 0;
        int idxMin  = 0;
        double maxVal=0.0;
        double minVal=10000.0;
        cv::Point ptSharp,ptSharp2,ptHead,ptHead2,ptTail;

        for (int j=0; j<X.size(); j++) {
           dXY[j] = (X[j]*X[j] + Y[j]*Y[j]);
           //if (dXY[j] > maxVal)
           //{
//               idxMax = j;
               maxVal = dXY[j];
  //         }

           if (dXY[j] < minVal) //Detect Tail
           {
               idxMin2 = idxMin;//Save As 2nd Smallest
               idxMin = j;
               minVal = dXY[j];
           }
        }

        ptSharp = curve[idxMin]; //Most Likely Tail Point
        ptSharp2 = curve[idxMin2]; //This Could Be Head

        ResampleCurve(smoothx,smoothy,resampledcurveX,resampledcurveY, gcFishContourSize,false);
        PolyLineMerge(curve,resampledcurveX,resampledcurveY);


        //Find Where Tail Point Is In the Resampled (Reduced) Contour
        int idxTail = findIndexClosesttoPoint(curve,ptSharp);

        //Put Tail Back to Curve In CAse it Got Smoothed Out
        std::vector<cv::Point>::iterator it = curve.begin();
        it += idxTail;
        curve.insert(it,ptSharp);

        //Copy with Tail At zero index
        vector<cv::Point> vcurveR( curve.begin() + idxTail ,curve.end());
        vcurveR.insert(vcurveR.end(),curve.begin(),curve.begin() + std::max(0,idxTail) );

        /// Find Head-Tail Point As Point Of Maximum Arc Lenght //
        //1st Approach is the noddy Way, n!
        int maxLen = 0;
        int lenA,lenB; //Arc Lenght ClockWise, And AntiClockWise
        int idxA,idxB,idxT;
        for (uint j=0; j < 2 ; j++) //Isolate to 1st found Tail Point
        {
            for (uint k=j+2; k< vcurveR.size(); k++)
            {
                vector<cv::Point>::const_iterator first = vcurveR.begin() + j;
                vector<cv::Point>::const_iterator last = vcurveR.begin() + k;
                vector<cv::Point>::const_iterator end = vcurveR.end();
                vector<cv::Point> vArc(first, last); //Copy Points Over And Test For arc Length
                vector<cv::Point> vArcRev(last, end); //Copy Points Over And Test For arc Length

                 lenA = cv::arcLength(vArc,false);
                 lenB = cv::arcLength(vArcRev,false); //Count Remaining Arc

                if (lenA >= lenB) //This Should Continuously Rise /
                 {
                    idxT = k;
                    break; //Crossed Over Mid Arc / So Distal point Found
                 }

            }//Search Across Points Ahead Of Starting Point

            //Is this the Longest Arc Found So Far? Save it , and Test new Starting Point
            if (lenA > maxLen)
            {
                maxLen = lenA;
                idxA = j; //Save the Curve Point Pair
                idxB = idxT;
            }

        }
        //
        //Check Which curve point is closest to the TailInflection Candidate - Initially Mark As Tail
        if (norm(curve[idxTail]-vcurveR[idxA]) < norm(curve[idxTail]-vcurveR[idxB]) )
        {
            ptTail = vcurveR[idxA];
            ptHead = vcurveR[idxB];
            ptHead2 = vcurveR[idxB-1];
        }
        else
        {
            ptTail = vcurveR[idxB];
            ptHead = vcurveR[idxA];
            ptHead2 = vcurveR[idxA-1];
        }
         //Verify Head Point is closest to the Centroid (COM) than tail / Otherwise Switch
        if (norm(ptTail-centroid)  < norm(ptHead-centroid))
        {
            cv::Point temp = ptTail;
            ptTail = ptHead;
            ptHead = temp;
        }


        //Replace Old Curve
        curve = vcurveR;

        /// \todo move this to some object prop, use ArcLength for fish Size Estimates
        gptHead = ptHead; //Hack To Get position For Template

        ///////////// END SMOOTHING

        ///\todo Make Contour Fish Like - Extend Tail ///
        //Find Tail Point- As the one with the sharpest Angle

        outfishbodycontours.push_back(curve);
        /////COMBINE - DRAW CONTOURS
        //Could Check if fishblob are contained (Doesn't matter if they are updated or not -
        // they should still fall within contour - )
        //cv::drawContours( maskFGImg, fgMaskcontours, kk, CV_RGB(0,0,0), cv::FILLED); //Erase Previous Fish Blob
        //Draw New One
        cv::drawContours( outFishMask, outfishbodycontours, (int)outfishbodycontours.size()-1, CV_RGB(255,255,255), cv::FILLED);
        cv::drawContours( outFishMask, outfishbodycontours, (int)outfishbodycontours.size()-1, CV_RGB(255,255,255),2); //Draw Thick Outline To Cover for Contour Losses
        cv::circle(outFishMask, (ptTail-ptHead)/6+ptTail,6,CV_RGB(255,255,255),cv::FILLED); //Add Trailing Expansion to the mask- In Case End bit of tail is not showing
        cv::circle(outFishMask, ptHead,4,CV_RGB(255,255,255),cv::FILLED);
         cv::circle(outFishMask, ptHead2,4,CV_RGB(255,255,255),cv::FILLED);

        //Erase Fish From Food Mask
        cv::drawContours( outFoodMask, outfishbodycontours, (int)outfishbodycontours.size()-1, CV_RGB(0,0,0),15);
        //cv::drawContours( outFoodMask, outfishbodycontours, (int)outfishbodycontours.size()-1, CV_RGB(0,0,0),4);

} //For Each Fish Contour




    //Merge Smoothed Contour Thresholded with BGMAsk //Add the masks so as to enhance fish features
    //cv::bitwise_or(outFishMask,maskFGImg,maskFGImg);
    //cv::bitwise_xor(outFishMask,outFoodMask,outFoodMask);
    //maskfishOnly.copyTo(maskFGImg);

    //threshold_output.copyTo(frameDebugD);

    if (bshowMask)
    {
        cv::imshow("Threshold Out",threshold_output);
        cv::imshow("Fish Mask",outFishMask);
        cv::imshow("Food Mask",outFoodMask); //Hollow Blobs For Detecting Food
        if (!fgMask.empty())
            cv::imshow("BG Model",fgMask);


    }

    // Release Should is done automatically anyway
    threshold_output.release();
    threshold_output_COMB.release();
   // maskfishOnly.release();
    //threshold_output_H.release();
}


///
/// \brief detectZfishFeatures - Used to create geometric representations of main zebrafish Features : Eyes, Body, tail
/// these are saved as point arrays on which angles and other measurements can be obtained
/// \param maskedGrayImg - IMage Masked so only fish is being shown Showing
/// \return
///
/// // \todo Optimize by re using fish contours already obtained in enhance fish mask
void detectZfishFeatures(MainWindow& window_main,const cv::Mat& fullImgIn,cv::Mat& fullImgOut,cv::Mat& imgFishHeadSeg, cv::Mat& maskedfishImg_gray, std::vector<std::vector<cv::Point> >& contours_body,std::vector<cv::Vec4i>& hierarchy_body)
{


//////No Longer Used Vars
    //        cv::RNG rng(12345);
    //    cv::Mat maskfishFeature,framelapl,framelapl_buffer;
//    cv::Mat frameCanny;
    //    cv::Mat grad,grad_x, grad_y;
    //    std::vector<std::vector<cv::Point> >hull( contours_body.size() );
    //std::vector<cv::Vec4i> hierarchy_canny; //Contour Relationships  [Next, Previous, First_Child, Parent]
    //std::vector<cv::Vec4i> hierarchy_laplace; //Contour Relationships  [Next, Previous, First_Child, Parent]
    /// Memory Crash on vector Init
    //std::vector<std::vector<cv::Point> > contours_laplace_clear; //For contours without markers
    //std::vector<cv::Vec4i> hierarchy_laplace_clear; //Contour Relationships  [Next, Previous, First_Child, Parent]
    //std::vector<std::vector<cv::Point> > fishfeatureContours( contours_laplace.size() );

/////////////////

    cv::Mat maskedImg_gray;
    cv::Mat maskedfishFeature_blur;
    // Memory Crash When Clearing Stack Here //
    //cv::Mat imgFishHeadSeg; //Thresholded / Or Edge Image Used In Detect Ellipses
    cv::Mat Mrot;

    //For Head Img//
    cv::Mat  imgFishAnterior,imgFishAnterior_Norm,imgFishHead,imgFishHeadProcessed; //imgTmp imgFishHeadEdge


    cv::Mat fullImg_colour;
    fullImgIn.convertTo(fullImg_colour,CV_8UC3);

    //fullImg_colour.copyTo(frameDebugC);

    /// Convert image to gray and blur it
    if (fullImgIn.depth() != CV_8U)
    {
        cv::cvtColor( fullImgIn, maskedImg_gray, cv::COLOR_BGR2GRAY );

    }
    else
        maskedImg_gray = fullImgIn; //Tautology

    cv::Mat fishTailFixed;

 ///Do not Use MaskedFish For Spine maskedfishImg_gray / + Fixed Contrast
    if (bUseMaskedFishForSpineDetect)
        fishTailFixed = maskedfishImg_gray*2.2;
    else
        fishTailFixed = maskedImg_gray*2.2;

    cv::GaussianBlur(fishTailFixed,maskedfishFeature_blur,cv::Size(11,11),7,7);

    //cv::imshow("BlugTail",maskedfishFeature_blur);
    //cv::imshow("ContrastTail",fishTailFixed);
    //Make image having masked all fish
    //maskedImg_gray.copyTo(maskedfishImg_gray,maskfishFGImg); //Mask The Laplacian //Input Already Masked


    ////Template Matching Is already Done On Fish Blob/Object
    //Pick The largest dimension and Make A Square
    cv::Size szTempIcon(std::max(gLastfishimg_template.cols,gLastfishimg_template.rows),std::max(gLastfishimg_template.cols,gLastfishimg_template.rows));
   // cv::Point rotCentre = cv::Point(szTempIcon.width/2,szTempIcon.height/2);


//    ///Detect Head Feature //
//    std::cout << "Match template on #fish:" << vfishmodels.size() << std::endl;
    for (fishModels::iterator it=vfishmodels.begin(); it!=vfishmodels.end(); ++it)
    //for (int z=0;z<vfishmodels.size();z++) //  fishModel* fish = vfishmodels[z];
    {

          fishModel* fish = (*it).second;

          //fish->bearingAngle   = AngleIdx;
            if (fish == 0 || fish->inactiveFrames > 1)
                continue;

          //Draw A general Region Where the FIsh Is located,
          cv::Point centre = fish->ptRotCentre; //top_left + rotCentre;
          //cv::Point centroid = fish->ptRotCentre ; // cv::Point2f(fish->track->centroid.x,fish->track->centroid.y);
          cv::Point pBound1 = cv::Point(max(0,min(maskedImg_gray.cols,centre.x-gFishBoundBoxSize)), max(0,min(maskedImg_gray.rows,centre.y-gFishBoundBoxSize)));
          cv::Point pBound2 = cv::Point(max(0,min(maskedImg_gray.cols,centre.x+gFishBoundBoxSize)), max(0,min(maskedImg_gray.rows,centre.y+gFishBoundBoxSize)));

          cv::Rect rectFish(pBound1,pBound2);

          //cv::rectangle(fullImgOut,rectFish,CV_RGB(20,200,150),2); //Identify Fish Region Bound In Cyan Square
#ifdef _ZTFDEBUG_

#endif
          // cv::Mat fishRegion(maskedImg_gray,rectFish); //Get Sub Region Image

          //0 Degrees Is along the Y Axis Looking Upwards

          ///Update Template Box Bound
          int bestAngleinDeg = fish->bearingAngle;
          cv::RotatedRect fishRotAnteriorBox(centre, cv::Size(gLastfishimg_template.cols,gLastfishimg_template.rows),bestAngleinDeg);
          /// Save Anterior Bound
          fish->bodyRotBound = fishRotAnteriorBox;

          //Locate Eyes In A box
          double lengthLine = 9;
          cv::Point2f ptEyeMid;

          ///Mark Point Between Eyes
          //Convert From Degrees and adjust to y Axis at 0 degrees (Ie flip of x,y)
          ptEyeMid.x =centre.x+lengthLine*sin((bestAngleinDeg)*(M_PI/180.0));
          ptEyeMid.y =centre.y-lengthLine*cos((bestAngleinDeg)*(M_PI/180.0)); //y=0 is the top left corner
          fish->midEyePoint = ptEyeMid;

          //Display MidEye Point
          //cv:circle(frameDebugC,ptEyeMid,1,CV_RGB(155,155,15),1);

          //Make A rectangle that surrounds part of the image that has been template matched
          cv::RotatedRect fishEyeBox(ptEyeMid, cv::Size(gLastfishimg_template.cols/2+3,gLastfishimg_template.cols/2+3),bestAngleinDeg);

          // Get Image Region Where the template Match occured
          //- Expand image so as to be able to fit the template When Rotated Orthonormally
          //Custom Bounding Box Needs to allow for RotRect To be rotated Orthonormally
          cv::Rect rectfishAnteriorBound = rectFish; //Use A square // fishRotAnteriorBox.boundingRect();
          cv::Size szFishAnteriorNorm(min(rectfishAnteriorBound.width,rectfishAnteriorBound.height)+4,max(rectfishAnteriorBound.width,rectfishAnteriorBound.height)+4); //Size Of Norm Image
          //Rot Centre Relative To Bounding Box Of UnNormed Image
          //cv::Point2f ptFishAnteriorRotCentre = (cv::Point2f)fishRotAnteriorBox.center-(cv::Point2f)rectfishAnteriorBound.tl();

          //Define Regions and Sizes for extracting Orthonormal Fish
          //Top Left Corner of templateSized Rect relative to Rectangle Centered in Normed Img
          cv::Size szTemplateImg = gLastfishimg_template.size();
          //cv::Point ptTopLeftTemplate(szFishAnteriorNorm.width/2-szTemplateImg.width/2,szFishAnteriorNorm.height/2-szTemplateImg.height/2);
          cv::Point ptTopLeftTemplate(rectfishAnteriorBound.width/2-szTemplateImg.width/2,rectfishAnteriorBound.height/2-szTemplateImg.height/2);


          cv::Rect rectFishTemplateBound = cv::Rect(ptTopLeftTemplate,szTemplateImg);
          cv::Size szHeadImg(min(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height),max(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height)*0.75);
//          cv::Point ptTopLeftHead(ptTopLeftTemplate.x,0);//(szFishAnteriorNorm.width/2-szTemplateImg.width/2,szFishAnteriorNorm.height/2-szTemplateImg.height/2);
          cv::Rect rectFishHeadBound = cv::Rect(ptTopLeftTemplate,szHeadImg);


          ///Make Normalized Fish View
           tEllipsoids vell;
           vell.clear();


           //maskedImg_gray.copyTo(imgTmp); //imgTmp Contain full frame Image in Gray
           //Threshold The Match Check Bounds Within Image
           cv::Rect imgBounds(0,0,maskedImg_gray.cols,maskedImg_gray.rows);

           if (!( //Looks Like a fish is found, now Check Bounds // gmaxVal > gTemplateMatchThreshold &&
               imgBounds.contains(rectfishAnteriorBound.br()) &&
                   imgBounds.contains(rectfishAnteriorBound.tl())))
               continue; //This Fish Is out Of Bounds /

              maskedImg_gray(rectfishAnteriorBound).copyTo(imgFishAnterior);
//              if (bUseEllipseEdgeFittingMethod)
//                frameCanny(rectfishAnteriorBound).copyTo(imgFishHeadEdge);
              //get Rotated Box Centre Coords relative to the cut-out of the anterior Body - This we use to rotate the image
              ///\note The centre of the Bounding Box could also do


              //cv::Point ptRotCenter = cv::Point(szFishAnteriorNorm.width/2,szFishAnteriorNorm.height/2);
              //cv::Point ptRotCenter = cv::Point(imgFishAnterior.cols/2,imgFishAnterior.rows/2);
              ///Make Rotation MAtrix cv::Point(imgFishAnterior.cols/2,imgFishAnterior.rows/2)
              cv::Point2f ptRotCenter = fishRotAnteriorBox.center - (cv::Point2f)rectfishAnteriorBound.tl();
             // ptRotCenter.x = ptRotCenter.x*cos(bestAngleinDeg*M_PI/180.0);
             // ptRotCenter.y = ptRotCenter.y*sin(bestAngleinDeg*M_PI/180.0);

              Mrot = cv::getRotationMatrix2D( ptRotCenter,bestAngleinDeg,1.0); //Rotate Upwards
              //cv::Mat Mrot = cv::getRotationMatrix2D(-fishRotHeadBox.center,bestAngleinDeg,1.0); //Rotate Upwards

              ///Make Rotation Transformation
              //Need to fix size of Upright/Normed Image
              cv::warpAffine(imgFishAnterior,imgFishAnterior_Norm,Mrot,szFishAnteriorNorm);

//if (bUseEllipseEdgeFittingMethod)
//              cv::warpAffine(imgFishHeadEdge,imgFishHeadEdge,Mrot,szFishAnteriorNorm);


              /// Store Norm Image as Template - If Flag Is set
              if (bStoreThisTemplate)
              {   std::stringstream ssMsg;
                  //Cut Down To Template Size
                  imgFishAnterior       = imgFishAnterior_Norm(rectFishTemplateBound);
                  addTemplateToCache(imgFishAnterior,gFishTemplateCache,gnumberOfTemplatesInCache);
                  //Try This New Template On the Next Search
                  iLastKnownGoodTemplateRow = gnumberOfTemplatesInCache-1;
                  fish->idxTemplateRow = iLastKnownGoodTemplateRow;
                  window_main.saveTemplateImage(imgFishAnterior);
                  ssMsg << "New Template Added, count is now #"<<gnumberOfTemplatesInCache << " NewRowIdx: " << iLastKnownGoodTemplateRow;
                  window_main.LogEvent(QString::fromStdString(ssMsg.str() ));
                  bStoreThisTemplate = false;
              }


              //Allow For Sub Optimal Matches To be processed Up to Here //
              //if (fish->templateScore < gTemplateMatchThreshold)
              //    continue; //Skip This Model Fish And Check the next one

              /// Prepare Norm Head Pic for Eye Detection Draw Centers for Reference and cleaner Masks
              //Draw  Rotation Centre of Transformation to Norm
              cv::circle(imgFishAnterior,ptRotCenter,3,CV_RGB(100,140,140),1);
              //cv::imshow("IsolatedAnterior",imgFishAnterior);

              ///Draw * Isolation Mask *  Of Eyes Around RotationCentre
              cv::Point ptMask(ptRotCenter.x,ptRotCenter.y+4);
//              cv::circle(imgFishAnterior_Norm,ptMask,giHeadIsolationMaskVOffset,CV_RGB(0,0,0),-1); //Mask Body
//              cv::line(imgFishAnterior_Norm,ptMask,cv::Point(imgFishAnterior_Norm.cols/2,0),CV_RGB(0,0,0),1);//Split Eyes
              imgFishHead           = imgFishAnterior_Norm(rectFishHeadBound);



              /// EYE DETECTION Report Results to Output Frame //
              int ret = detectEllipses(imgFishHead,vell,imgFishHeadSeg,imgFishHeadProcessed);
              std::stringstream ss;

            if (ret < 2)
            {
                ss << " Eye Detection Error - Check Threshold";
                window_main.LogEvent(QString::fromStdString(ss.str()));
                //fish->leftEyeTheta = 180;
                //fish->rightEyeTheta = 180;
                //fish->leftEye.fitscore = fish->rightEye.fitscore = 0;
                fish->nFailedEyeDetectionCount++;
                //std::clog << ss.str() << std::endl;
            }


              ///Paste Eye Processed Head IMage to Into Top Right corner of Larger Image
              cv::Rect rpasteregion(fullImgOut.cols-imgFishHeadProcessed.cols,0,imgFishHeadProcessed.cols,imgFishHeadProcessed.rows );
              //  show_histogram("HeadHist",imgFishHead);
              if (imgFishHeadProcessed.u)
                imgFishHeadProcessed.copyTo(fullImgOut(rpasteregion));

              //headImgOut = imgFishHeadSeg.clone(); //Return As INdividual Image Too which is then Shown On GUI Graphics Object
              //imgFishHeadSeg.copy
              //imgFishHeadSeg.release(); //Decrement Ref Counter


              /// Set Detected Eyes Back to Fish Features
              ///  Print Out Values
              /// \todo Figure out Why/how is it that nan Values Appeared in Output File : NA Values in ./Tracked07-12-17/LiveFed/Empty//AutoSet420fps_07-12-17_WTLiveFed4Empty_286_005_tracks_2.csv
              ss.str(""); //Empty String
              ss.precision(3);
              if (vell.size() > 0)
              {//Left Eye Detected First
                  tDetectedEllipsoid lEye = vell.at(0); //L Eye Is pushed 1st
                  fish->leftEye           = lEye;
                  fish->leftEyeTheta      = lEye.rectEllipse.angle;
                  ss << "L:" << fish->leftEyeTheta;
                  cv::putText(fullImgOut,ss.str(),cv::Point(rpasteregion.br().x-75,rpasteregion.br().y+10),CV_FONT_NORMAL,0.4,CV_RGB(250,250,0),1 );
              }else
              { //Set To Not detected
                  ss << "L Eye Detection Error - Check Threshold";
                  window_main.LogEvent(QString::fromStdString(ss.str()));

                  fish->leftEye       = tDetectedEllipsoid(cv::RotatedRect(),0);
                  fish->leftEyeTheta  = 180;
                  fish->leftEye.fitscore = 0;
                  fish->nFailedEyeDetectionCount++;
              }


              ss.str(""); //Empty String
              if (vell.size() > 1)
              {
                tDetectedEllipsoid rEye = vell.at(1); //R Eye Is pushed 2nd
                fish->rightEye          = rEye;
                fish->rightEyeTheta     = rEye.rectEllipse.angle;
                ss << "R:"  << fish->rightEyeTheta;
                cv::putText(fullImgOut,ss.str(),cv::Point(rpasteregion.br().x-75,rpasteregion.br().y+25),CV_FONT_NORMAL,0.4,CV_RGB(250,250,0),1 );
              }else
              { //Set To Not detected
                  ss << "R Eye Detection Error - Check Threshold";
                  window_main.LogEvent(QString::fromStdString(ss.str()));

                  fish->rightEye       = tDetectedEllipsoid(cv::RotatedRect(),0);
                  fish->rightEyeTheta  = 180;
                  fish->rightEye.fitscore = 0;
                  fish->nFailedEyeDetectionCount++;
              }


              ///If Both Eyes Detected Then Print Vergence Angle
              if (fish->leftEye.fitscore > 20 && fish->rightEye.fitscore > 20)
              {
                  ss.str(""); //Empty String
                  ss << "V:"  << fish->leftEyeTheta - fish->rightEyeTheta;
                  cv::putText(fullImgOut,ss.str(),cv::Point(rpasteregion.br().x-75,rpasteregion.br().y+40),CV_FONT_NORMAL,0.4,CV_RGB(250,250,0),1 );
                  fish->nFailedEyeDetectionCount = 0; //Reset Error Count
              }


              //Check If Too Many Eye Detection Failures - Then Switch Template
              if (fish->nFailedEyeDetectionCount > 10)
              {
                    fish->idxTemplateRow = iLastKnownGoodTemplateRow = (rand() % static_cast<int>(gnumberOfTemplatesInCache - 0 + 1));//Start From RANDOM rOW On Next Search
                    pwindow_main->LogEvent(QString("[warning] Too Many Eye detection Failures - Change Template Randomly to :" + QString::number(iLastKnownGoodTemplateRow)));
              }

              /// SPINE Fitting And Drawing ///
              /// \note two methods
              if (contours_body.size() > 0 && bFitSpineToTail)
              {
              //Look for Top Level Contour

               //Fit Spine to Countour MEthod
#ifdef _USEFITSPINETOCONTOUR
               int idxFish = findMatchingContour(contours_body,hierarchy_body,centre,2);
               if (idxFish>=0)
                fish->fitSpineToContour(maskedImg_gray,contours_body,0,idxFish);
#endif

               /// Main Method Uses Pixel Intensity //
               fish->fitSpineToIntensity(maskedfishFeature_blur,gFitTailIntensityScanAngleDeg);
               fish->drawSpine(fullImgOut);

               /// Optionally Check For Errors Using Spine to Contour Fitting Periodically//
#ifdef _USEPERIODICSPINETOCONTOUR_TEST
               if ( (pwindow_main->nFrame % (uint)gfVidfps) == 0)
               {
                   int idxFish = findMatchingContour(contours_body,hierarchy_body,centre,2);
                   if (idxFish>=0)
                        fish->fitSpineToContour(maskedImg_gray,contours_body,0,idxFish);

                   qDebug() << "Spine Tail Fit Error :" << fish->lastTailFitError;
               }

               //If Convergece TimedOut Then likely the fit is stuck with High Residual and no gradient
               //Best To reset Spine and Start Over Next Time
               /// \todo Make this parameter threshold formal
               if (abs(fish->lastTailFitError) > fish->c_fitErrorPerContourPoint)
               {
                   fish->resetSpine(); //No Solution Found So Reset
                   pwindow_main->LogEvent(QString("[warning] Reset Spine. lastTailFitError ") + QString::number(fish->lastTailFitError) + QString(" > c_fitErrorPerContourPoint") );

               }

#endif


                //cv::imshow("BlurredFish",maskedfishFeature_blur);
              }
             /// END OF Fit Spine ////



              //Eye Detection Ret > 0

    } //For eAch Fish Model


        bEyesDetected = false; //Flip Back to off in case it was eye features were marked for saving



    //Draw On Canny Img

#ifdef _ZTFDEBUG
    ///DEBUG show all contours -Edge
     frameCanny.convertTo(frameCanny, CV_8UC3);
    for( size_t i = 0; i< contours_canny.size(); i++ )
    {
         cv::drawContours( frameDebugC, contours_canny, (int)i, CV_RGB(200,0,60), 1,8,hierarchy_canny);
    }
    qDebug() << "maskedfishFeature_blur.u->refcount ==" << maskedfishFeature_blur.u->refcount;
    qDebug() << "Mrot.u->refcount ==" << Mrot.u->refcount;
    qDebug() << "imgFishAnterior.u->refcount ==" << imgFishAnterior.u->refcount;
    qDebug() << "imgFishAnterior_Norm.u->refcount ==" << imgFishAnterior_Norm.u->refcount;
    qDebug() << "imgFishHead.u->refcount == " << imgFishHead.u->refcount;
    qDebug() << "maskedImg_gray.u->refcount=" << maskedImg_gray.u->refcount;
    if (imgFishHeadProcessed.u)
        qDebug() << "imgFishHeadProcessed.u->refcount ==" << imgFishHeadProcessed.u->refcount;

#endif


    //cv::imshow("Edges Canny",frameCanny);
    //cv::imshow("Edges Laplace",framelapl);

    /// Clearn Up Check For Leaks //
    /// Local allocated cv:Mat must have 1 ref at this point



    //assert(maskedfishFeature_blur.u->refcount == 1);

    maskedfishFeature_blur.release();


    //assert(Mrot.u->refcount == 1);

    Mrot.release();


    //assert(imgFishAnterior.u->refcount == 1);

    imgFishAnterior.release();



    //assert(imgFishAnterior_Norm.u->refcount == 1);


    imgFishAnterior_Norm.release();


    //assert(imgFishHead.u->refcount == 1);

    imgFishHead.release();


    //assert(imgFishHeadProcessed.u->refcount == 1);
    imgFishHeadProcessed.release();

    //assert(maskedImg_gray.u->refcount == 2); //1 Ref Comes From InpuTIMgs


    maskedImg_gray.release();



} //DetectZFeatures


/**
* @function thresh_callback
*/
void thresh_callback(int, void* )
{

    if (g_BGthresh % 2 == 0)
        g_BGthresh ++;

    if (g_Segthresh <= 3) g_Segthresh = 3;

    if (g_Segthresh%2 == 0)
        g_Segthresh++;

    if (gi_CannyThres <2)
        gi_CannyThres = 2;

  //  Aperture size should be odd between 3 and 7 in function Canny
    if (gi_CannyThresSmall % 2 == 0)
        gi_CannyThresSmall ++;
    if (gi_CannyThresSmall <3)
        gi_CannyThresSmall =3;
    if (gi_CannyThresSmall >7)
        gi_CannyThresSmall =7;


    for (fishModels::iterator ft  = vfishmodels.begin(); ft!=vfishmodels.end(); ++ft)
    {
         fishModel* pfish = ft->second;
         pfish->c_spineSegL           = gFishTailSpineSegmentLength;


    }

//    if (!pGHT.empty())
//    {
//        pGHT->setCannyHighThresh(gi_CannyThres);
//        pGHT->setCannyLowThresh(gi_CannyThresSmall);

//        pGHTGuil->setScaleThresh(gi_VotesSThres);
//        pGHTGuil->setAngleThresh(gi_VotesAThres);
//        pGHTGuil->setPosThresh(gi_VotesPThres);

//        //Ptr<GeneralizedHoughBallard> ballard = static_cast<Ptr<GeneralizedHoughBallard>> pGHT;
//        //if (gi_ThresholdMatching>0)
//        //    pGHTBallard->setVotesThreshold(gi_ThresholdMatching);

//    }





}



//////////////////////////////////////////////////////////////////////////////
///
///
/// process_mem_usage(double &, double &) - takes two doubles by reference,
/// attempts to read the system-dependent data for a process' virtual memory
/// size and resident set size, and return the results in KB.
///
/// On failure, returns 0.0, 0.0

void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}


///////////////////// For Image Sequence /////////////


//unsigned int trackImageSequencefiles(MainWindow& window_main)
//{

//    cv::Mat frame,frameMasked,fgMask,outframe,outframeHead;

//    QString inVideoDirname = QFileDialog::getExistingDirectory(&window_main,"Select folder with video images to track", gstroutDirCSV);

//    //unsigned int istartFrame = 0;
//    unsigned int nFrame = 0;

//    QStringList strImageNames; //Save Passed Files Here

//    qDebug() << "Open File Sequence in : " << inVideoDirname;

//        ///* Obtain BG Model LOOP//////////////
//        //QDirIterator itBGimgDir(inVideoDirname, QDirIterator::Subdirectories);
//        QStringList fileFilters; fileFilters << "*.png" << "*.tiff" << "*.pgm" << "*.png";
//        QStringList imgFiles = QDir(inVideoDirname).entryList(fileFilters,QDir::Files,QDir::Name);
//        inVideoDirname.append('/');
//        QListIterator<QString> itfile (imgFiles);
//        while (itfile.hasNext() && !bExiting)
//        {
//          QString filename = itfile.next();
//          qDebug() << filename;
//          std::string filepath = filename.prepend(inVideoDirname ).toStdString();
//          //std::cout << filepath << std::endl;

//          frame  = cv::imread(filepath , CV_LOAD_IMAGE_COLOR);
//          //Contrast Brightness
//          //frame.convertTo(frame, -1, 1, 0); //increase the contrast (double)


//          nFrame++;
//          window_main.nFrame = nFrame;
//          if (!updateBGFrame(frame,fgMask,nFrame)) //Stop when BG learning says so
//            break;

//          /// Display Output //
//          frameMasked = cv::Mat::zeros(frameMasked.rows, frameMasked.cols, CV_8U);
//          frame.copyTo(frameMasked,fgMask);

//          ///Display Output
//          cv::imshow(gstrwinName,frameMasked);
//          //cv::displayOverlay(gstrwinName,"Press 'e' when features Have been detected" , 10000 );

//          window_main.showVideoFrame(frame,nFrame); //Show On QT Window
//          //cv::imshow(gstrwinName + " FG Mask", fgMask);
//          //Check For input Control
//          //keyboard = cv::waitKey( cFrameDelayms );
//          checkPauseRun(&window_main,keyboard,nFrame);
//        }


//    ///\brief LOOP Tracking Process images with Obtained BG Model - Now Start over images afresh
//    nFrame = 0;
//    window_main.nFrame = nFrame;

//    //Show Video list to process
//    //Go through Each Video - Hold Last Frame N , make it the start of the next vid.
//    std::cout << "Starting Tracking  processing" << std::endl;

//    itfile.toFront();
//    while (itfile.hasNext() && !bExiting)
//    {
//      QString filename = itfile.next();
//      qDebug() << filename;
//      std::string filepath = filename.prepend(inVideoDirname).toStdString();

//      //std::cout << filepath << std::endl;

//       frame  = cv::imread(filepath , CV_LOAD_IMAGE_COLOR);
//       if (!frame.data)
//       {
//            std::cerr << "Could not open next Image frame." << std::endl;
//            std::exit(EXIT_FAILURE);
//       }
//       //if (frame.depth() < 3) //Need To increase depth if we are to paint on this frame
//       //     cv::cvtColor( frame, frame, cv::COLOR_GRAY2RGB);

//       //Contrast Brightness
//       //frame.convertTo(frame, -1, 1.2, 0); //increase the contrast (double)
//       nFrame++;
//       window_main.nFrame = nFrame;

//              //Make Global Roi on 1st frame
//       if (nFrame == 1)
//       {
//           //Add Global Roi
//           ltROI newROI(cv::Point(frame.cols/2,frame.rows/2),ptROI2);
//           addROI(newROI);
//       }

//      // std::cout << " Now Processing : "<< itimgDir.fileName().toStdString() ;

//       /// Frame The Fish ///
//       frameMasked = cv::Mat::zeros(frame.rows, frame.cols,frame.type());
//       frame.copyTo(frameMasked,fgMask);


//       processFrame(window_main,frame,fgMask,nFrame,outframe,outframeHead);
//       cv::imshow(gstrwinName + " FishOnly",frameMasked);

//       /// Display Output //
//       window_main.showVideoFrame(outframe,nFrame); //Show On QT Window

//       if (bshowMask)
//       {
//            cv::imshow(gstrwinName + " FG Mask", fgMask);
//            //cv::imshow(gstrwinName + " FG Fish Mask", fgMaskFish); //The Circle mask
//       }


//       window_main.setWindowTitle("Tracking:" + filename);
//       //keyboard = cv::waitKey( cFrameDelayms );
//       checkPauseRun(&window_main,keyboard,nFrame);
//    }
//    return nFrame;
//}///trackImageSequencefiles



