#include <config.h>

#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include "cereal/types/vector.hpp"
#include <fstream>

// Gaussian Curve Smoothing Kernels For fish Contour//
//std::vector<double> gGaussian,dgGaussian,d2gGaussian;



QElapsedTimer gTimer;
QFile outfishdatafile;
QFile outfooddatafile;
QFile EyeDetectorRL; //Reinforcement Learned Behaviour For Eye Segmentation -

std::ofstream foutLog;//Used for Logging To File

//QString outfilename;
//std::string gstrwinName = "FishFrame";


//Global Matrices Used to show debug images
cv::Mat frameDebugA,frameDebugB,frameDebugC,frameDebugD;
cv::Mat gframeCurrent,gframeLast; //Updated in processVideo Global Var Holding Copy of current and previous frame - usefull for opticflows
cv::Mat gframeBGImage;

//Morphological Kernels
cv::Mat kernelOpen;
cv::Mat kernelDilateMOGMask;
cv::Mat kernelOpenfish;
cv::Mat kernelClose;

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
cv::Point gptHead,gptTail; //Candidate Fish Contour Position Of HEad - Use for template Detect

//For Inset Pasting
//cv::Rect rect_pasteregion;

//uint gi_MaxFishID;
//uint gi_MaxFoodID; //Declared in Model Header Files

class MainWindow;
extern MainWindow pwindow_main;
extern trackerState gTrackerState;



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
    pMOG2 =  cv::createBackgroundSubtractorMOG2(gTrackerState.MOGhistory, 20,false);
#endif

    pMOG2->setHistory(gTrackerState.MOGhistory);
    pMOG2->setNMixtures(20);
    pMOG2->setBackgroundRatio(gTrackerState.gdMOGBGRatio);

}
/// end of bg substractor init //


/// \brief Load internal and external template images memory cache //
int initDetectionTemplates()
{
    cv::Mat lastfish_template_img;// OUr Fish Image Template
    ///////////////////////////////////////
    /// Setup Fish Body Template Cache //
    int idxTempl;

    for (idxTempl=0; idxTempl<gTrackerState.nTemplatesToLoad;idxTempl++)
    {
        trackerState::loadFromQrc(QString::fromStdString(gTrackerState.strTemplateImg + to_string(idxTempl+1) + std::string(".pgm")),lastfish_template_img,IMREAD_GRAYSCALE); //  loadImage(strTemplateImg);
        if (lastfish_template_img.empty())
        {
            std::cerr << "Could not load template" << std::endl;
            exit(-1);
        }
        //Add to Global List Of Template Images
        gTrackerState.vTemplImg.push_back(lastfish_template_img);
        //Add to Cache and generate all All Angle Varations
        addTemplateToCache(lastfish_template_img,gFishTemplateCache,idxTempl); //Increments Index
    }

    // Set Template Size
    gszTemplateImg.width = lastfish_template_img.size().width; //Save TO Global Size Variable
    gszTemplateImg.height = lastfish_template_img.size().height; //Save TO Global Size Variable

    // Set Paster Region for Inset Image
    gTrackerState.rect_pasteregion.x = (640-gszTemplateImg.width*2);
    gTrackerState.rect_pasteregion.y = 0;
    gTrackerState.rect_pasteregion.width = gszTemplateImg.width*2; //For the upsampled image
    gTrackerState.rect_pasteregion.height = gszTemplateImg.height*2;

    gTrackerState.gstroutDirTemplates = gTrackerState.gstroutDirCSV + ("/templates/");
    int ifileCount = loadTemplatesFromDirectory(QString::fromStdString(  gTrackerState.gstroutDirTemplates) );

    //Make Mean Fish And Add to Cache
     cv::Mat templFrame = makeMeanTemplateImage(gTrackerState.vTemplImg);
     addTemplateToCache(templFrame,gFishTemplateCache,gTrackerState.gnumberOfTemplatesInCache);
     gTrackerState.gLastfishimg_template = templFrame; //Set To Global
 #if defined(_DEBUG)
     cv::imshow("Template Cache",gFishTemplateCache);
#endif

    return (gTrackerState.gnumberOfTemplatesInCache);
    /// END OF FISH TEMPLATES ///
}



/// State Class Methods //

trackerState::trackerState()
{
    bPaused     = true;
    bROIChanged = true;
    bPaused     = false;
    bshowMask   = false;
    bTracking   = true; //Start By Tracking by default
    bExiting    = false;



}

void trackerState::setVidFps(float fps)
{
    gfVidfps                      = fps;
    gcMaxFoodModelInactiveFrames  = gfVidfps*2; //Number of frames inactive (Not Matched to a Blob) until track is deleted
    gcMinFoodModelActiveFrames    = gfVidfps/5;
    gFoodReportInterval           = gfVidfps; //Report Food every second

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



/// \brief Process user provided config params and set corresponding internal/global variables
void trackerState::initGlobalParams(cv::CommandLineParser& parser,QStringList& inVidFileNames)
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
        gstroutDirCSV  = QString::fromStdString(soutFolder).toStdString();

    }
    else
    {
      outfilename  = QFileDialog::getSaveFileName(nullptr, "Save tracks to output","VX_pos.csv", "CSV files (*.csv);", nullptr, nullptr); // getting the filename (full path)
      gstroutDirCSV = outfilename.left(outfilename.lastIndexOf("/")).toStdString();
    }


    std::cout << "Csv Output Dir is " << gstroutDirCSV  << "\n " <<std::endl;


 /// Check if vid file provided in arguments.
 /// If File exists added to video file list,
 /// otherwise save directory and open dialogue to choose a file from there
    if (parser.has("invideofile"))
    {   QString fvidFileName = QString::fromStdString( parser.get<std::string>("invideofile") );
        QFileInfo ovidfile(fvidFileName ) ;

        if ( ovidfile.absoluteDir().exists()) //Check if vid file exists before appending to list
            gstrinDirVid = ovidfile.absoluteDir().absolutePath().toStdString();
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
    kernelOpen          = cv::getStructuringElement(cv::MORPH_CROSS,cv::Size(1,1),cv::Point(-1,-1));
    kernelDilateMOGMask = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(3,3),cv::Point(-1,-1));
    kernelOpenfish      = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(3,3),cv::Point(-1,-1)); //Note When Using Grad Morp / and Low res images this needs to be 3,3
    kernelClose         = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(5,5),cv::Point(-1,-1));

    /// create Gaussian Smoothing kernels for Contour //
    assert(gTrackerState.dGaussContourKernelSize % 2 == 1); //M is an odd number
    getGaussianDerivs(dGaussContourKernelSigma,dGaussContourKernelSize,gGaussian,dgGaussian,d2gGaussian);

}
/// END OF INIT GLOBAL PARAMS //

/// \brief Initializes ROI at start of tracking depending on user params / either large circle or user defined/configurable polygon
void  trackerState::initROI()
{
    //Rect Roi Keep Away from L-R Edges to Avoid Tracking IR lightRing Edges
    ptROI1 = cv::Point(gFishBoundBoxSize*2+1,gFishBoundBoxSize/2);
    ptROI2 = cv::Point(640-gFishBoundBoxSize*2,gFishBoundBoxSize/2);
    ptROI3 = cv::Point(640-gFishBoundBoxSize*2,512-gFishBoundBoxSize/2);
    ptROI4 = cv::Point(gFishBoundBoxSize*2+1,512-gFishBoundBoxSize/2);


    /// Init Polygon ROI ///
    ///Make A Rectangular Roi Default //
    if ( bMakeCustomROIRegion)
    {
        std::vector<cv::Point> vPolygon;
        vPolygon.push_back(ptROI1); vPolygon.push_back(ptROI2); vPolygon.push_back(ptROI3); vPolygon.push_back(ptROI4);
        ltROI rectRoi(vPolygon);
        vRoi.push_back(rectRoi);
     }
      else //Make Default ROI Region
    {
        ptROI2.x = (640- gFishBoundBoxSize)/2;
        ptROI2.y = gszTemplateImg.height/3;
    //Add Global Roi - Center - Radius
        ltROI newROI(cv::Point(640/2,520/2),ptROI2);
        addROI(newROI);
    }
}

/// \brief Load Q Resources
void trackerState::loadFromQrc(QString qrc,cv::Mat& imRes,int flag )
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
