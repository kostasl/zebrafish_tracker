#include <config.h>


#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include "cereal/types/vector.hpp"
#include <fstream>

/// VIDEO AND BACKGROUND PROCESSING //
float gfVidfps                  = 410;
const unsigned int MOGhistory   = 100;//Use 100 frames Distributed across the video length To Find What the BGModel is
double gdMOGBGRatio             = 0.05; ///If a foreground pixel keeps semi-constant value for about backgroundRatio*history frames, it's considered background and added to the model as a center of a new component.

//Processing Loop delay
uint cFrameDelayms              = 1;

const double dLearningRate                = 1.0/(MOGhistory); //Learning Rate During Initial BG Modelling - Learn Slow So 1st Playbacl Frame doesnt look new anymore
const double dLearningRateNominal         = 0.05; //Fast Rate as BG Learning Allows for threshold after BGSubstract operation to Augment The mask
double dBGMaskAccumulateSpeed             = 1.0/(4.0*MOGhistory);

/// BLOB DETECTION Filters //
//Area Filters
const double dMeanBlobArea                    = 100; //Initial Value that will get updated
const double dVarBlobArea                     = 20;
const unsigned int gc_fishLength            = 100; //px Length Of Fish
const unsigned int thresh_fishblobarea      = 350; //Min area above which to Filter The fish blobs
const unsigned int thresh_maxfishblobarea = 2550; //Min area above which to Filter The fish blobs
const unsigned int gthres_maxfoodblobarea = thresh_fishblobarea/3;


/// Constants ///
int gcMaxFishModelInactiveFrames    = 150; //Number of frames inactive until track is deleted
int gcMaxFoodModelInactiveFrames    = gfVidfps*2; //Number of frames inactive (Not Matched to a Blob) until track is deleted
int gcMinFoodModelActiveFrames      = gfVidfps/10; //Min Number of consecutive frames it needs to be active  otherwise its deleted
const int gMaxClusterRadiusFoodToBlob   = 6;
const int thActive                      = 0;// Deprecated If a track becomes inactive but it has been active less than thActive frames, the track will be deleted.
const int gc_FishTailSpineSegmentLength_init = 9;
//const int thDistanceFish                = 250; //NOT USED / Blobs need to overlap Threshold for distance between track (last known point) -to blob assignement
//const int thDistanceFood                = 4; //Threshold for distance between track-to blob assignement
int gFoodReportInterval           = (int)gfVidfps;

const int nTemplatesToLoad  = 19; //Number of Templates To Load Into Cache - These need to exist as images in QtResources


///Segmentation Params
int g_Segthresh             = 15; //Applied On THe BG Substracted Image / Image Threshold to segment BG - Fish Segmentation uses a higher 2Xg_Segthresh threshold
int g_SegFoodThesMin        = 16; //Low thres For Food Detection / Doing Gradual Step Wise with SimpleBlob
int g_SegFoodThesMax        = g_Segthresh+5; //Up thres Scan For Food Detection / Doing Gradual Step Wise with SimpleBlob
int g_SegInnerthreshMult    = 3; //Image Threshold for Inner FIsh Features //Deprecated
int g_BGthresh              = 10; //BG threshold segmentation
int gi_ThresholdMatching    = 10; /// Minimum Score to accept that a contour has been found
int gi_FoodModelNumberLimit = 250; /// Maximum Number of Food Objects /Prey To track
bool gOptimizeShapeMatching = false; ///Set to false To disable matchShapes in FindMatching Contour

int gi_CannyThres           = 150;
int gi_CannyThresSmall      = 50; //Aperture size should be odd between 3 and 7 in function Canny
int gi_maxEllipseMajor      = 21; /// thres  for Eye Ellipse Detection methods
int gi_minEllipseMajor      = 10; ///thres for Eye Ellipse Detection methods (These Values Tested Worked Best)
int gi_minEllipseMinor      = 0; /// ellipse detection width - When 0 it allows for detecting straight line
int giEyeIsolationMaskRadius = 10; ///Mask circle between eyes
int gi_VotesEllipseThres        = 9; //Votes thres for The Backup Ellipse Detection Based on the Hough Transform
int gthresEyeSeg                = 3; //Additional Adjustment for Adaptive Threshold  For Eye Segmentation In Isolated Head IMage
int gthresEyeSegL               = 2;
int gnumberOfTemplatesInCache   = 0; //INcreases As new Are Added
float gDisplacementThreshold    = 2.0; //Distance That Fish Is displaced so as to consider active and Record A point For the rendered Track /
int gFishBoundBoxSize           = 24; /// pixel width/radius of bounding Box When Isolating the fish's head From the image
int gFishTailSpineSegmentLength     = gc_FishTailSpineSegmentLength_init;



int gFitTailIntensityScanAngleDeg   = 40; //
const int gFishTailSpineSegmentCount= ZTF_TAILSPINECOUNT;
const int gcFishContourSize         = ZTF_FISHCONTOURSIZE;
const int gMaxFitIterations         = ZTF_TAILFITMAXITERATIONS; //Constant For Max Iteration to Fit Tail Spine to Fish Contour

int giHeadIsolationMaskVOffset      = 24; //Vertical Distance to draw  Mask and Threshold Sampling Arc in Fish Head Mask

///Fish Features Detection Params
int gFishTemplateAngleSteps     = 1;
int gEyeTemplateAngleSteps      = 5;
int iEyeMaskSepWidth            = 25; //5 px width vertical line separates the eyes for segmentation
double eyeStepIncrement         = 0.1;
double gTemplateMatchThreshold  = 0.89; //If not higher than 0.9 The fish body can be matched at extremeties
int iLastKnownGoodTemplateRow   = 0;
int iFishAngleOffset            = 0;
double gUserReward              = 0; //User feedback for reinforcement learning
int iTemplateMatchFailCounter   = 0; //Counts the number of consecutive times template failed to match
//using namespace std;


// Global Control Vars ///

int keyboard; //input from keyboard
int screenx,screeny;
bool bshowMask; //True will show the BGSubstracted IMage/Processed Mask
bool bROIChanged;
bool bPaused;
bool bStartPaused;
bool bExiting;
bool bTracking;
bool bTrackFood         = true;
bool bAddPreyManually   = false;
bool bMeasure2pDistance = true; /// A mode allowing 2point distance measurement
bool bTrackFish         = true;
bool bRecordToFile      = true;
bool bSaveImages            = false;
bool b1stPointSet;
bool bMouseLButtonDown;
//bool bSaveBlobsToFile; //Check in fnct processBlobs - saves output CSV
bool bEyesDetected = false; ///Flip True to save eye shape feature for future detection
bool bStoreThisTemplate             = false;
bool bDraggingTemplateCentre        = false;
bool bUseEllipseEdgeFittingMethod   = true; //Allow to Use the 2nd Efficient Method of Ellipsoid Fitting if the 1st one fails - Set to false to Make trakcing Faster
bool bFitSpineToTail                = true; // Runs The Contour And Tail Fitting Spine Optimization Algorith
bool bStartFrameChanged         = false; /// When True, the Video Processing loop stops /and reloads video starting from new Start Position

bool bRenderToDisplay           = true; ///Updates Screen to User When True
bool bRenderWithAlpha           = false;
bool bDrawFoodBlob              = false; ///Draw circle around identified food blobs (prior to model matching)
bool bOffLineTracking           = false; ///Skip Frequent Display Updates So as to  Speed Up Tracking
bool bBlindSourceTracking       = false; /// Used for Data Labelling, so as to hide the data source/group/condition
bool bStaticAccumulatedBGMaskRemove       = false; /// Remove Pixs from FG mask that have been shown static in the Accumulated Mask after the BGLearning Phase
bool bUseBGModelling                      = true; ///Use BG Modelling TO Segment FG Objects
bool gbUpdateBGModel                      = true; //When Set a new BGModel Is learned at the beginning of the next video
bool gbUpdateBGModelOnAllVids             = true; //When Set a new BGModel Is learned at the beginning of the next video
bool bApplyFishMaskBeforeFeatureDetection = true; ///Pass the masked image of the fish to the feature detector
bool bSkipExisting                        = false; /// If A Tracker DataFile Exists Then Skip This Video
bool bMakeCustomROIRegion                 = false; /// Uses Point array to construct
bool bUseMaskedFishForSpineDetect         = true; /// When True, The Spine Is fit to the Masked Fish Image- Which Could Be problematic if The contour is not detected Well
bool bTemplateSearchThroughRows           = false; /// Stops TemplateFind to Scan Through All Rows (diff temaplte images)- speeding up search + fail - Rows still Randomly Switch between attempts
bool bRemovePixelNoise                    = false; //Run Gaussian Filter Noise Reduction During Tracking
bool bUseGPU                              = false;
bool bUseOpenCL                           = true;
bool bUseHistEqualization                 = true; //To enhance to contrast in Eye Ellipse detection
/// \todo Make this path relative or embed resource
//string strTemplateImg = "/home/kostasl/workspace/cam_preycapture/src/zebraprey_track/img/fishbody_tmp.pgm";
std::string strTemplateImg = ":/img/fishbody_tmp"; ///Load From Resource

uint uiStartFrame = 1;
uint uiStopFrame = 0;


// Gaussian Curve Smoothing Kernels For fish Contour//
std::vector<double> gGaussian,dgGaussian,d2gGaussian;
double sigma = 6.0;
int M = static_cast<int>(round((18.0*sigma+1.0) / 2.0)) * 2 - 1; //Gaussian Kernel Size




// Other fonts:
//   CV_FONT_HERSHEY_SIMPLEX, CV_FONT_HERSHEY_PLAIN,
//   CV_FONT_HERSHEY_DUPLEX, CV_FONT_HERSHEY_COMPLEX,
//   CV_FONT_HERSHEY_TRIPLEX, CV_FONT_HERSHEY_COMPLEX_SMALL,
//   CV_FONT_HERSHEY_SCRIPT_SIMPLEX, CV_FONT_HERSHEY_SCRIPT_COMPLEX
int trackFnt = CV_FONT_HERSHEY_SIMPLEX;  //Font for Reporting - Tracking
float trackFntScale = 0.6f;


QElapsedTimer gTimer;
QFile outfishdatafile;
QFile outfooddatafile;
QFile EyeDetectorRL; //Reinforcement Learned Behaviour For Eye Segmentation -

std::ofstream foutLog;//Used for Logging To File

QString outfilename;
std::string gstrwinName = "FishFrame";
QString gstroutDirCSV,gstrinDirVid,gstrvidFilename; //The Output Directory



//Global Matrices Used to show debug images
cv::Mat frameDebugA,frameDebugB,frameDebugC,frameDebugD;
cv::Mat gframeCurrent,gframeLast; //Updated in processVideo Global Var Holding Copy of current and previous frame - usefull for opticflows
cv::Mat gframeBGImage;

//Morphological Kernels
cv::Mat kernelOpen;
cv::Mat kernelDilateMOGMask;
cv::Mat kernelOpenfish;
cv::Mat kernelClose;
cv::Mat gLastfishimg_template;// OUr Fish Image Template
cv::Mat gFishTemplateCache; //A mosaic image contaning copies of template across different angles
//cv::Mat gEyeTemplateCache; //A mosaic image contaning copies of template across different angles


//Global CUda Utility Matrices Used to Tranfser Images To GPU
// Defined here to save reallocation Time
#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
        cv::cuda::GpuMat dframe_mask; //Passed to MOG Cuda
        cv::cuda::GpuMat dframe_gray; // For Denoising
        cv::cuda::GpuMat dframe_thres; // Used In Mask Enhancement
        cv::Ptr<cv::cuda::TemplateMatching> gpu_MatchAlg;// For Template Matching
        Ptr<cuda::Filter> gpu_DilateFilter;
#endif



cv::Size gszTemplateImg;

//cv::Ptr<cv::BackgroundSubtractor> pMOG; //MOG Background subtractor
//cv::Ptr<cv::BackgroundSubtractorKNN> pKNN; //MOG Background subtractor
//cv::Ptr<cv::bgsegm::BackgroundSubtractorGMG> pGMG; //GMG Background subtractor

// Fish Detection //
cv::Ptr<cv::BackgroundSubtractorMOG2> pMOG2; //MOG2 Background subtractor
//cv::Ptr<cv::GeneralizedHough> pGHT;
//cv::Ptr<cv::GeneralizedHoughBallard> pGHTBallard;
//cv::Ptr<cv::GeneralizedHoughGuil> pGHTGuil;

/// \todo using a global var is a quick hack to transfer info from blob/Mask processing to fishmodel / Need to change the Blob Struct to do this properly
cv::Point gptHead; //Candidate Fish Contour Position Of HEad - Use for template Detect

ltROIlist vRoi;
//
//Rect Roi Keep Away from L-R Edges to Avoid Tracking IR lightRing Edges
cv::Point ptROI1 = cv::Point(gFishBoundBoxSize*2+1,gFishBoundBoxSize/2);
cv::Point ptROI2 = cv::Point(640-gFishBoundBoxSize*2,gFishBoundBoxSize/2);
cv::Point ptROI3 = cv::Point(640-gFishBoundBoxSize*2,512-gFishBoundBoxSize/2);
cv::Point ptROI4 = cv::Point(gFishBoundBoxSize*2+1,512-gFishBoundBoxSize/2);

//For Inset Pasting
cv::Rect rect_pasteregion;

uint gi_MaxFishID;
uint gi_MaxFoodID; //Declared in Model Header Files

class MainWindow;
extern MainWindow pwindow_main;




/// \brief Load Q Resources
void loadFromQrc(QString qrc,cv::Mat& imRes,int flag = cv::IMREAD_COLOR)
{
    //double tic = double(getTickCount());

    QFile file(qrc);

    if(file.open(QIODevice::ReadOnly))
    {
        qint64 sz = file.size();
        std::vector<uchar> buf(sz);
        file.read((char*)buf.data(), sz);
        imRes = cv::imdecode(buf, flag);
    }else
        std::cerr << " Could not load template file " << qrc.toStdString();

    //double toc = (double(getTickCount()) - tic) * 1000.0 / getTickFrequency();
    //qDebug() << "OpenCV loading time: " << toc;

}

/// MAIN FUNCTION - ENTRY POINT ////

/// \brief Count Number of different Characters Between str1 and str2
int compString(QString str1,QString str2)
{
    int ret =0;
  for (int j=0;j<std::min(str1.length(),str2.length());j++)
        if (str1.mid(j,1) != str2.mid(j,1))
            ret++;

  return ret;
}


/// \brief sigsev Error handlers picking up unhandled errors
void installErrorHandlers()
{


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
       fprintf(stderr, "Error setting signal handler for %d (%s)\n",
         SIGSEGV, strsignal(SIGSEGV));

       exit(EXIT_FAILURE);
      }
 /// ERROR HANDLER SIGSEV /////////
}

/// \brief Process user provided config params and set corresponding internal/global variables
void initGlobalParams(cv::CommandLineParser& parser,QStringList& inVidFileNames)
{

    QString outfilename;

    if (parser.has("help"))
    {
        parser.printMessage();
        return;
    }

    if (parser.has("outputdir"))
    {
        std::string soutFolder   = parser.get<std::string>("outputdir");
        std::clog << "Cmd Line OutDir : " << soutFolder <<std::endl;
        gstroutDirCSV  = QString::fromStdString(soutFolder);

    }
    else
    {
      outfilename  = QFileDialog::getSaveFileName(nullptr, "Save tracks to output","VX_pos.csv", "CSV files (*.csv);", nullptr, nullptr); // getting the filename (full path)
      gstroutDirCSV = outfilename.left(outfilename.lastIndexOf("/"));
    }


    std::cout << "Csv Output Dir is " << gstroutDirCSV.toStdString()  << "\n " <<std::endl;


 /// Check if vid file provided in arguments.
 /// If File exists added to video file list,
 /// otherwise save directory and open dialogue to choose a file from there
    if (parser.has("invideofile"))
    {   QString fvidFileName = QString::fromStdString( parser.get<std::string>("invideofile") );
        QFileInfo ovidfile(fvidFileName ) ;

        if ( ovidfile.absoluteDir().exists()) //Check if vid file exists before appending to list
            gstrinDirVid = ovidfile.absoluteDir().absolutePath();
        else
            gstrinDirVid = gstroutDirCSV; //Set Def. Dir for dialogue to the outDir

        if (ovidfile.exists() && ovidfile.isFile())
            inVidFileNames.append( fvidFileName );
    }


    if (parser.has("invideolist"))
    {
        qDebug() << "Load Video File List " <<  QString::fromStdString(parser.get<std::string>("invideolist"));
        QFile fvidfile( QString::fromStdString(parser.get<std::string>("invideolist")) );
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
        qDebug() << "Set Log File To " <<  QString::fromStdString( parser.get<std::string>("logtofile") );

        QFileInfo oLogPath( QString::fromStdString(parser.get<std::string>("logtofile") ) );
        if (!oLogPath.absoluteDir().exists())
            QDir().mkpath(oLogPath.absoluteDir().absolutePath()); //Make Path To Logs

        foutLog.open(oLogPath.absoluteFilePath().toStdString());
         // Set the rdbuf of clog.
         std::clog.rdbuf(foutLog.rdbuf());
         std::cerr.rdbuf(foutLog.rdbuf());
    }

    // Read In Flag To enable Fish Tracking / FishBlob Processing
    if (parser.has("TrackFish"))
        bTrackFish = (parser.get<int>("TrackFish") == 1)?true:false;

    //Check If We Are BG Modelling / BEst to switch off when Labelling Hunting Events
    if (parser.has("ModelBG"))
        bUseBGModelling = (parser.get<int>("ModelBG") == 1)?true:false;

    if (parser.has("ModelBGOnAllVids"))
        gbUpdateBGModelOnAllVids = (parser.get<int>("ModelBGOnAllVids") == 1)?true:false;

    if (parser.has("SkipTracked"))
        bSkipExisting = (parser.get<int>("SkipTracked") == 1)?true:false;

    if (parser.has("PolygonROI"))
        bMakeCustomROIRegion = (parser.get<int>("PolygonROI") == 1)?true:false;

      if (parser.has("BGThreshold"))
          g_Segthresh = parser.get<int>("BGThreshold");

    if (parser.has("FilterPixelNoise"))
    {
        bRemovePixelNoise = (parser.get<int>("FilterPixelNoise") == 1)?true:false;
        std::clog << "Remove Pixel Noise Filter Is On" << std::endl;
    }

    if (parser.has("startpaused"))
        bStartPaused = (parser.get<int>("startpaused") == 1)?true:false;

    if (parser.has("HideDataSource"))
        bBlindSourceTracking = (parser.get<int>("HideDataSource") == 1)?true:false;

    if (parser.has("EyeHistEqualization"))
        bUseHistEqualization = (parser.get<int>("EyeHistEqualization") == 1)?true:false;

    ///Disable OPENCL in case SEG Fault is hit - usually from MOG when running multiple tracker processes
    if (parser.has("DisableOpenCL")){
        if (parser.get<int>("DisableOpenCL") == 1)
        {
            cv::ocl::setUseOpenCL(false);
            bUseOpenCL =false;
        }else{
            cv::ocl::setUseOpenCL(true);
            bUseOpenCL =true;
        }
    }

    if (parser.has("EnableCUDA"))
           bUseGPU = (parser.get<int>("EnableCUDA") == 1)?true:false;

    if (parser.has("MeasureMode"))
        bMeasure2pDistance = (parser.get<int>("MeasureMode") == 1)?true:false;


    uiStartFrame = parser.get<uint>("startframe");
    uiStopFrame = parser.get<uint>("stopframe");

    ///* Create Morphological Kernel Elements used in processFrame *///
    kernelOpen      = cv::getStructuringElement(cv::MORPH_CROSS,cv::Size(1,1),cv::Point(-1,-1));
    kernelDilateMOGMask = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(3,3),cv::Point(-1,-1));
    kernelOpenfish  = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(3,3),cv::Point(-1,-1)); //Note When Using Grad Morp / and Low res images this needs to be 3,3
    kernelClose     = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(5,5),cv::Point(-1,-1));

    /// create Gaussian Smoothing kernels for Contour //
    getGaussianDerivs(sigma,M,gGaussian,dgGaussian,d2gGaussian);

}
/// END OF INIT GLOBAL PARAMS //

/// \brief Initializes ROI at start of tracking depending on user params / either large circle or user defined/configurable polygon
void initROI()
{
    /// Init Polygon ROI ///
    ///Make A Rectangular Roi Default //
    if (bMakeCustomROIRegion)
    {
        std::vector<cv::Point> vPolygon;
        vPolygon.push_back(ptROI1); vPolygon.push_back(ptROI2); vPolygon.push_back(ptROI3); vPolygon.push_back(ptROI4);
        ltROI rectRoi(vPolygon);
        vRoi.push_back(rectRoi);
     }
      else //Make Default ROI Region
    {
        ptROI2.x = (640-gFishBoundBoxSize)/2;
        ptROI2.y = gszTemplateImg.height/3;
    //Add Global Roi - Center - Radius
        ltROI newROI(cv::Point(640/2,520/2),ptROI2);
        addROI(newROI);
    }
}

/// \brief Initialize BG substractor objects, depending on options / can use cuda
void initBGSubstraction()
{

    /// CUDA Version Of BG MOG //
#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    if (bUseGPU)
    {
        pMOG2 = cv::cuda::createBackgroundSubtractorMOG2(MOGhistory,20,false);
        gpu_MatchAlg = cv::cuda::createTemplateMatching(CV_8U, CV_TM_CCORR_NORMED);
        gpu_DilateFilter = cuda::createMorphologyFilter(MORPH_DILATE, CV_8U, kernelDilateMOGMask);
    }else
        pMOG2 =  cv::createBackgroundSubtractorMOG2(MOGhistory, 20,false);
#else
    //Doesn't matter if cuda FLAG is enabled
    pMOG2 =  cv::createBackgroundSubtractorMOG2(MOGhistory, 20,false);
#endif

    pMOG2->setHistory(MOGhistory);
    pMOG2->setNMixtures(20);
    pMOG2->setBackgroundRatio(gdMOGBGRatio); ///

}
/// end of bg substractor init //


/// \brief Load internal and external template images memory cache //
int initDetectionTemplates()
{
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

    // Set Template Size
    gszTemplateImg.width = gLastfishimg_template.size().width; //Save TO Global Size Variable
    gszTemplateImg.height = gLastfishimg_template.size().height; //Save TO Global Size Variable

    // Set Paster Region for Inset Image
    rect_pasteregion.x = (640-gszTemplateImg.width*2);
    rect_pasteregion.y = 0;
    rect_pasteregion.width = gszTemplateImg.width*2; //For the upsampled image
    rect_pasteregion.height = gszTemplateImg.height*2;

    int ifileCount = loadTemplatesFromDirectory(gstroutDirCSV + QString("/templates/"));
    return (ifileCount+nTemplatesToLoad);
    /// END OF FISH TEMPLATES ///
}

/// State Class Methods //

trackerState::trackerState()
{

}


void trackerState::saveState(std::string strFilename)
{
    std::ofstream os(strFilename);
    cereal::XMLOutputArchive archive(os);
    this->serialize(archive); //save State Values

    os.flush();

}

void trackerState::loadState(std::string strFilename)
{

    /// Load Archived values if they Exists
    /// Load Saved Learned Behaviour
     assert(strFilename > 0);
     qDebug() << "Load tracker State:" << QString::fromStdString(strFilename);
     std::ifstream is(strFilename);
     if (is.is_open())
     {

       try
         {
           cereal::XMLInputArchive archive(is);
           archive(userROI); //Load State Value

         }catch (QString e)
         {
                 qDebug() << "Failed to open Tracker State file:" << e;
         }


     }
}
