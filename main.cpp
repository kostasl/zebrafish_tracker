///*
//// 25/11/2015 : kostasl Testing OpenCV bg substraction - to process larva in vial recording timelapses.
 //// App uses BG substyraction MOG2, with a slow learning rate.
 //// then Uses Open and Close / Dilation contraction techniques to get rid of noise and fill gaps
 //// Then uses cvBlob library to track and count larva on BG-processed images.
 //// The lib sources have been included to the source and slightly modified in update tracks to fix a bug.
 ////
 ///* User:
 ///* Chooses input video file, then on the second dialogue choose the text file to export track info in CSV format.
 ///* The green box defines the region over which the larvae are counted-tracked and recorded to file.
 ///* Once the video begins to show, use to left mouse clicks to define a new region in the image over which you want to count the larvae.
 ///* Press p to pause Image. once paused:
 ///*  s to save snapshots in CSV outdir pics subfolder.
 ///*  2 Left Clicks to define the 2 points of region-of interest for tracking.
 ///*  m to show the masked image of the larva against BG.
 ///*  t Start Tracking
 ///*  q Exit Quit application
 ///*
 ///* NOTE: ChFanging ROI hits SEG. FAULTs in update tracks of the library. So I made setting of ROI only once.
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
 ////////


#include <larvatrack.h>
#include <ellipse_detect.h>
#include <template_detect.h>
#include <zfttracks.h>

#include <QDirIterator>
#include <QDir>
#include <QDebug>
//#include <QThread>
#include <QTime>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"

//#include <opencv2/bgsegm.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/video/background_segm.hpp>

#include <GUI/mainwindow.h>

#include <CSS/CurveCSS.h> ///Curve Smoothing and Matching

/// Constants ///
const int inactiveFrameCount            = 30000; //Number of frames inactive until track is deleted
const int thActive                      = 0;// If a track becomes inactive but it has been active less than thActive frames, the track will be deleted.
const int thDistanceFish                = 150; //Threshold for distance between track-to blob assignement
const int thDistanceFood                = 15; //Threshold for distance between track-to blob assignement
const double dLearningRateNominal       = 0.000;

/// Vars With Initial Values  -
//Area Filters
double dMeanBlobArea                    = 100; //Initial Value that will get updated
double dVarBlobArea                     = 20;
const unsigned int gc_fishLength        = 100; //px Length Of Fish
const unsigned int thresh_fishblobarea  = 220; //Min area above which to Filter The fish blobs

//BG History
float gfVidfps              = 300;
const int MOGhistory        = gfVidfps*2;
//Processing Loop delay
uint cFrameDelayms          = 1;
double dLearningRate        = 1.0/(2*MOGhistory);

///Segmentation Params
int g_Segthresh             = 35; //Image Threshold to segment the FIsh/ Body and Tail
int g_SegInnerthreshMult    = 3; //Image Threshold for Inner FIsh Features //Deprecated
int g_BGthresh              = 10; //BG threshold segmentation
int gi_ThresholdMatching    = 10; /// Minimum Score to accept that a contour has been found
bool gOptimizeShapeMatching = false; ///Set to false To disable matchShapes in FindMatching Contour
int gi_CannyThres           = 150;
int gi_CannyThresSmall      = 50; //Aperture size should be odd between 3 and 7 in function Canny
int gi_maxEllipseMajor      = 11; // thres for Hough Transform
int gi_minEllipseMajor      = 7; //thres for Hough Transform
int gi_VotesEllipseThres    = 9; //Votes thres for Hough Transform
int gthresEyeSeg            = 105; //Threshold For Eye Segmentation In Isolated Head IMage
int gnumberOfTemplatesInCache  = 0; //INcreases As new Are Added
const int nTemplatesToLoad      = 5; //Number of Templates To Load Into Cache - These need to exist as images in QtResources
float gDisplacementThreshold = 0.5; //Distance That Fish Is displaced so as to consider active and Record A point For the rendered Track /
int gFishBoundBoxSize        = 20; /// pixel width/radius of bounding Box When Isolating the fish's head From the image

///Fish Features Detection Params
int gFishTemplateAngleSteps     = 2;
int gEyeTemplateAngleSteps      = 5;
double gMatchShapeThreshold     = 0.65;
int iLastKnownGoodTemplateRow   = 0;
int iLastKnownGoodTemplateCol   = 0;
//using namespace std;


///Global Variables
QElapsedTimer gTimer;
QString outfilename;
std::string gstrwinName = "FishFrame";
QString gstroutDirCSV; //The Output Directory

//Global Matrices Used to show debug images
cv::Mat frameDebugA,frameDebugB,frameDebugC,frameDebugD;

//cv::Ptr<cv::BackgroundSubtractor> pMOG; //MOG Background subtractor
cv::Ptr<cv::BackgroundSubtractorMOG2> pMOG2; //MOG2 Background subtractor
//cv::Ptr<cv::BackgroundSubtractorKNN> pKNN; //MOG Background subtractor
//cv::Ptr<cv::bgsegm::BackgroundSubtractorGMG> pGMG; //GMG Background subtractor

//Fish Detection
Ptr<GeneralizedHough> pGHT;
Ptr<GeneralizedHoughBallard> pGHTBallard;
Ptr<GeneralizedHoughGuil> pGHTGuil;

//Morphological Kernels
cv::Mat kernelOpen;
cv::Mat kernelOpenLaplace;
cv::Mat kernelOpenfish;
cv::Mat kernelClose;
cv::Mat fishbodyimg_template;// OUr Fish Image Template
cv::Mat gFishTemplateCache; //A mosaic image contaning copies of template across different angles
cv::Mat gEyeTemplateCache; //A mosaic image contaning copies of template across different angles


//Global Shortcut of Type conversion to legacy IplImage
//IplImage framefishMaskImg;


ltROI Circle( cv::Point(0,0) , cv::Point(1024,768));
ltROIlist vRoi;
cv::Point ptROI1 = cv::Point(320,240);
cv::Point ptROI2 = cv::Point(1,131);


//Structures to hold blobs & Tracks
//Blobs as identified by BG Substractions
//cvb::CvBlobs blobs; //All Blobs - Updated Ids on everyframe done by cvLabel function
//cvb::CvBlobs fishblobs;
//cvb::CvBlobs foodblobs;
cvb::CvTracks fishtracks;
cvb::CvTracks foodtracks;
cvb::CvTracks tracks; ///All tracks

//The fish ones are then revaluated using simple thresholding to obtain more accurate contours
fishModels vfishmodels; //Vector containing live fish models


// Other fonts:
//   CV_FONT_HERSHEY_SIMPLEX, CV_FONT_HERSHEY_PLAIN,
//   CV_FONT_HERSHEY_DUPLEX, CV_FONT_HERSHEY_COMPLEX,
//   CV_FONT_HERSHEY_TRIPLEX, CV_FONT_HERSHEY_COMPLEX_SMALL,
//   CV_FONT_HERSHEY_SCRIPT_SIMPLEX, CV_FONT_HERSHEY_SCRIPT_COMPLEX
int trackFnt = cv::FONT_HERSHEY_COMPLEX_SMALL;  //Font for Reporting - Tracking
float trackFntScale = 0.7;

// Global Control Vars ///

int keyboard; //input from keyboard
int screenx,screeny;
bool bshowMask; //True will show the BGSubstracted IMage/Processed Mask
bool bROIChanged;
bool bPaused;
bool bExiting;
bool bTracking;
bool bSaveImages = false;
bool b1stPointSet;
bool bMouseLButtonDown;
bool bSaveBlobsToFile; //Check in fnct processBlobs - saves output CSV
bool bEyesDetected = false; ///Flip True to save eye shape feature for future detection
bool bStoreThisTemplate = false;
bool bDraggingTemplateCentre = false;

/// \todo Make this path relative or embed resource
//string strTemplateImg = "/home/kostasl/workspace/cam_preycapture/src/zebraprey_track/img/fishbody_tmp.pgm";
string strTemplateImg = ":/img/fishbody_tmp"; ///Load From Resource

static Mat loadImage(const string& name)
{
    Mat image = imread(name, IMREAD_GRAYSCALE);
    if (image.empty())
    {
        cerr << "Can't load image - " << name << endl;
        exit(-1);
    }
    return image;
}

Mat loadFromQrc(QString qrc, int flag = IMREAD_COLOR)
{
    //double tic = double(getTickCount());

    QFile file(qrc);
    Mat m;
    if(file.open(QIODevice::ReadOnly))
    {
        qint64 sz = file.size();
        std::vector<uchar> buf(sz);
        file.read((char*)buf.data(), sz);
        m = imdecode(buf, flag);
    }else
        std::cerr << " Could not load template file " << qrc.toStdString();

    //double toc = (double(getTickCount()) - tic) * 1000.0 / getTickFrequency();
    //qDebug() << "OpenCV loading time: " << toc;

    return m;
}

int main(int argc, char *argv[])
{
    bROIChanged = false;
    bPaused = false;
    bshowMask = false;
    bTracking = false;
    bExiting    = false;

    QApplication app(argc, argv);
    //QQmlApplicationEngine engine;
    MainWindow window_main;

    window_main.show();

    //window_main.showFullScreen();

    //outfilename.truncate(outfilename.lastIndexOf("."));
    QString outfilename = QFileDialog::getSaveFileName(0, "Save tracks to output","VX_pos.csv", "CSV files (*.csv);", 0, 0); // getting the filename (full path)
    //QString outDir = outfilename.left(outfilename.lastIndexOf('/') ).toStdString().c_str();
    gstroutDirCSV = outfilename.left(outfilename.lastIndexOf("/"));
    std::cout << "Csv Output Dir is " << gstroutDirCSV.toStdString()  << "\n " <<std::endl;

    // get the applications dir pah and expose it to QML
    //engine.load(QUrl(QStringLiteral("qrc:///main.qml")));

    gTimer.start();
    //create GUI windows


    //cv::namedWindow(gstrwinName,CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    //cv::namedWindow(gstrwinName + " FishOnly",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    cv::namedWindow("Debug A",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    cv::namedWindow("Debug B",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    cv::namedWindow("Debug C",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    cv::namedWindow("Debug D",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    //cv::namedWindow("Ellipse fit",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    cv::namedWindow("HeadHist",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);

    frameDebugA = cv::Mat::zeros(640, 480, CV_8U);
    frameDebugB = cv::Mat::zeros(640, 480, CV_8U);
    frameDebugC = cv::Mat::zeros(640, 480, CV_8U);
    frameDebugD = cv::Mat::zeros(640, 480, CV_8U);

    //set the callback function for any mouse event
    //cv::setMouseCallback(gstrwinName, CallBackFunc, NULL);

    cv::setMouseCallback("HeadHist", CallBackHistFunc, NULL);


    //Initialize The Track and blob vectors
    cvb::cvReleaseTracks(tracks); //Releasing All tracks will delete all track Objects
    //cvb::cvReleaseBlobs(blobs);
    ReleaseFishModels(vfishmodels);


    /// create Background Subtractor objects
    //(int history=500, double varThreshold=16, bool detectShadows=true
    //OPENCV 3

    pMOG2 =  cv::createBackgroundSubtractorMOG2(MOGhistory,16,false);
    pMOG2->setNMixtures(150);
    pMOG2->setBackgroundRatio(0.99);

    double dmog2TG = pMOG2->getVarThresholdGen();
    //pMOG2->setVarThresholdGen(1.0);
    double dmog2VT = pMOG2->getVarThreshold();
    pMOG2->setVarThreshold(3.0);

    ///////////////////////////////////////
    /// Setup Fish Body Template Cache //

    int idxTempl;

    for (idxTempl=0; idxTempl<nTemplatesToLoad;idxTempl++)
    {
        fishbodyimg_template = loadFromQrc(QString::fromStdString(strTemplateImg + to_string(idxTempl+1) + std::string(".pgm")),IMREAD_GRAYSCALE); //  loadImage(strTemplateImg);
        if (fishbodyimg_template.empty())
        {
            std::cerr << "Could not load template" << std::endl;
            exit(-1);
        }

        addTemplateToCache(fishbodyimg_template,gFishTemplateCache,idxTempl); //Increments Index
    }


    /// END OF FISH TEMPLATES ///

    ///Make TrackBars ///
    cv::createTrackbar( "Laplace Size:",  "Debug D", &g_BGthresh, 31.0, thresh_callback );
    cv::createTrackbar( "Fish Threshold:", "Debug D", &g_Segthresh, 151.0, thresh_callback );
    cv::createTrackbar( "Vote Threshold:", "Debug D", &gi_ThresholdMatching, 120.0, thresh_callback );
    cv::createTrackbar( "Canny Thres:", "Debug D", &gi_CannyThres, 350, thresh_callback );
    cv::createTrackbar( "Canny Thres Small:", "Debug D", &gi_CannyThresSmall, 100, thresh_callback );
    cv::createTrackbar( "Max Ellipse","Debug D", &gi_maxEllipseMajor, 20.0, thresh_callback );
    cv::createTrackbar( "Min Ellipse","Debug D", &gi_minEllipseMajor,10, thresh_callback );
    cv::createTrackbar( "Ellipse Votes:","Debug D", &gi_VotesEllipseThres, 800, thresh_callback );

    thresh_callback( 0, 0 );
    ///////////////

    //double mog2CThres = pMOG2->getComplexityReductionThreshold(); ///This parameter defines the number of samples needed to accept to prove the component exists. CT=0.05 is a default value for all the samples. By setting CT=0 you get an algorithm very similar to the standard Stauffer&Grimson algorithm.
    //pMOG2->setComplexityReductionThreshold(0.0);

    //(int history=200, int nmixtures=5, double backgroundRatio=0.7, double noiseSigma=0)
    //pMOG =   cv::bgsegm::createBackgroundSubtractorMOG(MOGhistory,12,0.05,0.00); //MOG approach
     //pKNN = cv::createBackgroundSubtractorKNN(MOGhistory,50,false);
//    pGMG =   cv::bgsegm::createBackgroundSubtractorGMG(MOGhistory,0.3); //GMG approach

    ///* Create Morphological Kernel Elements used in processFrame *///
    kernelOpen      = cv::getStructuringElement(cv::MORPH_CROSS,cv::Size(1,1),cv::Point(-1,-1));
    kernelOpenLaplace = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(1,1),cv::Point(-1,-1));
    kernelOpenfish  = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(3,3),cv::Point(-1,-1)); //Note When Using Grad Morp / and Low res images this needs to be 3,3
    kernelClose     = cv::getStructuringElement(cv::MORPH_ELLIPSE,cv::Size(1,1),cv::Point(-1,-1));


    //unsigned int hWnd = cvGetWindowHandle(sgstrwinName);
//    try{ //Check If cv is compiled with QT support //Remove otherwise
//        cv::setWindowTitle(strwinName, outfilename.toStdString());

//     //trackVideofiles(window_main);

//    }catch(int e)
//    {
//        std::cerr << "OpenCV not compiled with QT support! can display overlay" <<std::endl;
//    }

    //trackImageSequencefiles(window_main);

    trackVideofiles(window_main,outfilename);
    //destroy GUI windows
    cv::destroyAllWindows();
    //cv::waitKey(0);                                          // Wait for a keystroke in the window

    //pMOG->~BackgroundSubtractor();
    //pMOG2->~BackgroundSubtractor();
    //pKNN->~BackgroundSubtractor();
    //pGMG->~BackgroundSubtractor();

    //Empty The Track and blob vectors
    cvb::cvReleaseTracks(tracks);
    //cvb::cvReleaseBlobs(blobs);



    std::cout << "Total processing time : mins " << gTimer.elapsed()/60000.0 << std::endl;
///Clean Up //
    frameDebugA.deallocate();
    frameDebugB.deallocate();
    frameDebugC.deallocate();
    frameDebugD.deallocate();

    //app.quit();
    window_main.close();



    return app.exec();

}



unsigned int trackVideofiles(MainWindow& window_main,QString outputFile)
{
    cv::Mat fgMask;
    QString invideoname = "*.mpg";
    unsigned int istartFrame = 0;
    QStringList invideonames =QFileDialog::getOpenFileNames(0, "Select timelapse video to Process",gstroutDirCSV.toStdString().c_str(), "Video file (*.mpg *.avi *.mp4 *.h264 *.mkv *.tiff *.png *.jpg *.pgm)", 0, 0);

    //Show Video list to process
    std::cout << "Video List To process:" <<std::endl;
    for (int i = 0; i<invideonames.size() && i < 10; ++i)
    {
       invideoname = invideonames.at(i);
       std::cout << "*" <<  invideoname.toStdString() << std::endl;
    }

    //Go through Each Image/Video - Hold Last Frame N , make it the start of the next vid.
    for (int i = 0; i<invideonames.size(); ++i)
    {

       invideoname = invideonames.at(i);
       std::cout << " Now Processing : "<< invideoname.toStdString() <<std::endl;
       //cv::displayOverlay(gstrwinName,"file:" + invideoname.toStdString(), 10000 );

       //getBGModelFromVideo(fgMask, window_main,invideoname,outfilename,istartFrame);

       std::cout << "Press p to pause Video processing" << std::endl;

       istartFrame = processVideo(fgMask,window_main,invideoname,outputFile,istartFrame);

       window_main.setWindowTitle("Tracking:" + invideoname);

        if (istartFrame == 0)
        {
            std::cerr << "Could not load last video - Exiting loop." <<std::endl;
            break;
        }
    }
    return istartFrame;
}


unsigned int trackImageSequencefiles(MainWindow& window_main)
{

    cv::Mat frame,frameMasked,fgMask,outframe;
    QString inVideoDirname = QFileDialog::getExistingDirectory(&window_main,"Select folder with video images to track", gstroutDirCSV);

    unsigned int istartFrame = 0;
    unsigned int nFrame = 0;

    QStringList strImageNames; //Save Passed Files Here

    qDebug() << "Open File Sequence in : " << inVideoDirname;

        ///* Obtain BG Model LOOP//////////////
        //QDirIterator itBGimgDir(inVideoDirname, QDirIterator::Subdirectories);
        QStringList fileFilters; fileFilters << "*.png" << "*.tiff" << "*.pgm" << "*.png";
        QStringList imgFiles = QDir(inVideoDirname).entryList(fileFilters,QDir::Files,QDir::Name);
        inVideoDirname.append('/');
        QListIterator<QString> itfile (imgFiles);
        while (itfile.hasNext() && !bExiting)
        {
          QString filename = itfile.next();
          qDebug() << filename;
          std::string filepath = filename.prepend(inVideoDirname ).toStdString();
          //std::cout << filepath << std::endl;

          frame  = cv::imread(filepath , CV_LOAD_IMAGE_COLOR);
          //Contrast Brightness
          //frame.convertTo(frame, -1, 1, 0); //increase the contrast (double)


          nFrame++;
          window_main.nFrame = nFrame;
          if (!updateBGFrame(frame,fgMask,nFrame)) //Stop when BG learning says so
            break;

          /// Display Output //
          frameMasked = cv::Mat::zeros(frameMasked.rows, frameMasked.cols, CV_8U);
          frame.copyTo(frameMasked,fgMask);

          ///Display Output
          cv::imshow(gstrwinName,frameMasked);
          //cv::displayOverlay(gstrwinName,"Press 'e' when features Have been detected" , 10000 );

          window_main.showVideoFrame(frame,nFrame); //Show On QT Window
          //cv::imshow(gstrwinName + " FG Mask", fgMask);
          //Check For input Control
          //keyboard = cv::waitKey( cFrameDelayms );
          checkPauseRun(&window_main,keyboard,nFrame);
        }


    ///\brief LOOP Tracking Process images with Obtained BG Model - Now Start over images afresh
    nFrame = 0;
    window_main.nFrame = nFrame;

    //Show Video list to process
    //Go through Each Video - Hold Last Frame N , make it the start of the next vid.
    std::cout << "Starting Tracking  processing" << std::endl;

    itfile.toFront();
    while (itfile.hasNext() && !bExiting)
    {
      QString filename = itfile.next();
      qDebug() << filename;
      std::string filepath = filename.prepend(inVideoDirname).toStdString();

      //std::cout << filepath << std::endl;

       frame  = cv::imread(filepath , CV_LOAD_IMAGE_COLOR);
       if (!frame.data)
       {
            std::cerr << "Could not open next Image frame." << std::endl;
            std::exit(EXIT_FAILURE);
       }
       //if (frame.depth() < 3) //Need To increase depth if we are to paint on this frame
       //     cv::cvtColor( frame, frame, cv::COLOR_GRAY2RGB);

       //Contrast Brightness
       //frame.convertTo(frame, -1, 1.2, 0); //increase the contrast (double)
       nFrame++;
       window_main.nFrame = nFrame;

              //Make Global Roi on 1st frame
       if (nFrame == 1)
       {
           //Add Global Roi
           ltROI newROI(cv::Point(frame.cols/2,frame.rows/2),ptROI2);
           addROI(newROI);
       }

      // std::cout << " Now Processing : "<< itimgDir.fileName().toStdString() ;

       /// Frame The Fish ///
       frameMasked = cv::Mat::zeros(frame.rows, frame.cols,frame.type());
       frame.copyTo(frameMasked,fgMask);


       processFrame(frame,fgMask,frameMasked,nFrame,outframe);
       cv::imshow(gstrwinName + " FishOnly",frameMasked);

       /// Display Output //
       window_main.showVideoFrame(outframe,nFrame); //Show On QT Window

       if (bshowMask)
       {
            cv::imshow(gstrwinName + " FG Mask", fgMask);
            //cv::imshow(gstrwinName + " FG Fish Mask", fgMaskFish); //The Circle mask
       }


       window_main.setWindowTitle("Tracking:" + filename);
       //keyboard = cv::waitKey( cFrameDelayms );
       checkPauseRun(&window_main,keyboard,nFrame);
    }
    return nFrame;
}///trackImageSequencefiles

///*
///Create FG Model Image - Since target objects can be/will be are moving from the 1st frame, we need a statistical model
/// of the BG precalculated
///
unsigned int getBGModelFromVideo(cv::Mat& fgMask,MainWindow& window_main,QString videoFilename,QString outFileCSV,unsigned int startFrameCount)
{
        cv::Mat frame;
        unsigned int nFrame         = startFrameCount; //Current Frame Number

        std::cout << "Starting Background Model processing..." << std::endl;
        //create the capture object
        cv::VideoCapture capture(videoFilename.toStdString());
        if(!capture.isOpened())
        {
            //error in opening the video input
            std::cerr << "Unable to open video file: " << videoFilename.toStdString() << std::endl;
            std::exit(EXIT_FAILURE);
        }

        //read input data. ESC or 'q' for quitting
        while( !bExiting && (char)keyboard != 27 && nFrame <= MOGhistory)
        {
            //read the current frame
            if(!capture.read(frame))
            {
                if (nFrame == startFrameCount)
                {
                    std::cerr << "Unable to read first frame." << std::endl;
                    nFrame = 0; //Signals To caller that video could not be loaded.
                    exit(EXIT_FAILURE);
                }
                else
                {
                    std::cerr << "Unable to read next frame. So this video Is done." << std::endl;
                   std::cout << nFrame << " frames of Video processed. Move on to next timelapse video? " <<std::endl;
                  //  break;
                   continue;
               }
            }
            //Add frames from Last video
            nFrame = capture.get(CV_CAP_PROP_POS_FRAMES) + startFrameCount;
            window_main.nFrame = nFrame;
            window_main.tickProgress();

            /// Call Update BG Model ///
            updateBGFrame(frame,fgMask,nFrame);

            //Hold A copy of Frame With all txt
            //frame.copyTo(frameMasked);

            //cvb::CvBlobs blobs;
            //show the current frame and the fg masks
            //cv::imshow(gstrwinName, frame);
            //window_main.showVideoFrame(frame,nFrame); //Show On QT Window

            //cv::imshow(gstrwinName + " FG Mask", fgMask);
            //cv::imshow("FG Mask MOG", fgMaskMOG);
            //cv::imshow("FG Mask GMG ", fgMaskGMG);

           // if (!bTracking)
           //get the input from the keyboard
           //keyboard = cv::waitKey( cFrameDelayms );


           checkPauseRun(&window_main,keyboard,nFrame);


        } //main While loop
        //delete capture object
        capture.release();

        //delete kernel;
        //delete kernelClose;


        std::cout << "Background Processing  loop. Finished" << std::endl;

        return nFrame;
} ///trackImageSequencefile


void processFrame(cv::Mat& frame,cv::Mat& fgMask,cv::Mat& frameMasked, unsigned int nFrame,cv::Mat& outframe)
{
    std::vector<std::vector<cv::Point> > fishbodycontours;
    std::vector<cv::Vec4i> fishbodyhierarchy;
    IplImage lplframe;

    unsigned int nLarva         =  0;
    unsigned int nFood          =  0;
    double dblRatioPxChanged    =  0.0;

    std::string frameNumberString = to_string(nFrame);

    //For Morphological Filter
    ////cv::Size sz = cv::Size(3,3);
    //frame.copyTo(inputframe); //Keep Original Before Painting anything on it
    //update the background model
    //OPEN CV 2.4
    // dLearningRate is now Nominal value
    frame.copyTo(outframe); //Make Replicate On which we draw output
    //No Need For MOG!
    //pMOG2->apply(outframe, fgMask,dLearningRateNominal);


    ///DRAW ROI
    drawROI(outframe);


    //lplframe = frameMasked; //Convert to legacy format

    //cvb::CvBlobs blobs;
    ///DO Tracking
    if (bTracking)
    {
       //Simple Solution was to Use Contours To measure LUarvae

        //Draw THe fish Masks more accuratelly by threshold detection - Enhances full fish body detection
    //    enhanceFishMask(outframe, fgMask,fishbodycontours,fishbodyhierarchy);// Add fish Blobs
        cv::Mat fgFishMask,fgFishImgMasked;
        cv::Mat fgFoodMask,fgFoodImgMasked;
        enhanceMask(frame,fgMask,fgFishMask,fgFoodMask,fishbodycontours, fishbodyhierarchy);
        //frameMasked = cv::Mat::zeros(frame.rows, frame.cols,CV_8UC3);
        outframe.copyTo(fgFishImgMasked,fgFishMask); //Use Enhanced Mask
        //outframe.copyTo(fgFoodImgMasked,fgFoodMask); //Use Enhanced Mask
        //show the current frame and the fg masks
        //cv::imshow(gstrwinName + " FishOnly",frameMasked);



       // Filters Blobs between fish and food - save into global vectors
        //processBlobs(&lplframe,fgMask, blobs,tracks,gstroutDirCSV,frameNumberString,dMeanBlobArea);
        std::vector<cv::KeyPoint> ptFishblobs;
        processFishBlobs(fgFishImgMasked,fgFishMask, outframe , ptFishblobs);
        nLarva = ptFishblobs.size();


        cv::Mat maskedImg_gray,maskedfishImg_gray;
        /// Convert image to gray and blur it
        cv::cvtColor( frame, maskedImg_gray, cv::COLOR_BGR2GRAY );
        ////Make image having masked all fish
        //maskedImg_gray.copyTo(maskedfishImg_gray,fgMask); //Mask The Laplacian //Input Already Masked
        //Update Fish Models Against Image and Tracks
        //UpdateFishModels(maskedImg_gray,vfishmodels,fishtracks);


        UpdateFishModels(maskedImg_gray,vfishmodels,ptFishblobs,nFrame);
        //If A fish Is Detected Then Draw Its tracks
        fishModels::iterator ft = vfishmodels.begin();
        if (ft != vfishmodels.end())
        {
            fishModel* pfish = ft->second;
            assert(pfish);
            zftRenderTrack(pfish->zTrack, frame, outframe,CV_TRACK_RENDER_PATH , trackFnt );
        }

        ///Keep A Global List of all tracks?
        //Combine Lists into Tracks before rendering
//        tracks.clear();
        //tracks.insert(foodtracks.begin(),foodtracks.end() );
        //saveTracks(vfishmodels,trkoutFileCSV,frameNumberString);

        detectZfishFeatures(frame,outframe,fgMask,fishbodycontours,fishbodyhierarchy); //Creates & Updates Fish Models

        ///////  Process Food Blobs ////
        // Process Food blobs
        std::vector<cv::KeyPoint> ptFoodblobs;
        nFood = processFoodBlobs(fgFoodMask,fgFoodMask, outframe , ptFoodblobs); //Use Just The Mask




        if (bSaveImages)
        {
            saveImage(to_string(nFrame),gstroutDirCSV,outframe);
            cv::putText(frameDebugC, "Save ON", cv::Point(15, 600),
                    cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));

        }

    } //If Tracking

    fishbodycontours.clear();
    fishbodyhierarchy.clear();
    //Save to Disk

    ///

    ///TEXT INFO Put Info TextOn Frame
    //Frame Number
    std::stringstream ss;
    cv::rectangle(outframe, cv::Point(10, 2), cv::Point(100,20),
               CV_RGB(10,10,10), -1);
    cv::putText(outframe, frameNumberString,  cv::Point(15, 15),
            trackFnt, trackFntScale ,  CV_RGB(250,250,0));

    //Count on Original Frame
    std::stringstream strCount;
    strCount << "Nf:" << (nLarva) << " Nr:" << nFood;
    cv::rectangle(outframe, cv::Point(10, 25), cv::Point(80,45),  CV_RGB(10,10,10), -1);
    cv::putText(outframe, strCount.str(), cv::Point(15, 38),
           trackFnt, trackFntScale ,  CV_RGB(250,250,0));


    char buff[100];
    static double vm, rss;

    //Report Time
    std::sprintf(buff,"t: %0.2f",gTimer.elapsed()/(1000.0*60.0) );
    //strLearningRate << "dL:" << (double)(dLearningRate);
    cv::rectangle(outframe, cv::Point(10, 50), cv::Point(50,70), cv::Scalar(10,10,10), -1);
    cv::putText(outframe, buff, cv::Point(15, 63),
            trackFnt, trackFntScale , CV_RGB(250,250,0));

    //Time Rate - conv from ms to minutes
    ///Memory Usage
    if (nFrame%30)
    {
        //THats In KiB units /So 1Million is A Gigabyte
        process_mem_usage(vm, rss);
        //std::cout << "VM: " << vm/1024.0 << "; RSS: " << rss/1024.0 << endl;

    }
    std::sprintf(buff,"Vm: %0.2fMB;Rss:%0.2fMB",vm/1024.0,rss/1024.0);
    cv::rectangle(outframe, cv::Point(5, 490), cv::Point(80,510), cv::Scalar(10,10,10), -1);
    cv::putText(outframe, buff, cv::Point(10, 505),
            trackFnt,trackFntScale , CV_RGB(10,250,0));

} //End Of Process Frame


///
/// \brief updateBGFrame Update BG model for a fixed number of frames
/// \param frame
/// \param fgMask
/// \param nFrame
/// \return returns false when limit of updates is reached
///
bool updateBGFrame(cv::Mat& frame, cv::Mat& fgMask, unsigned int nFrame)
{

    bool ret = true;
    //Speed that stationary objects are removed
    double dblRatioPxChanged    = 0.0;

    //update the background model
    //OPEN CV 2.4
    if (nFrame > MOGhistory)
    {
        dLearningRate =dLearningRateNominal; //Nominal
        ret = false;
    }
    dblRatioPxChanged = (double)cv::countNonZero(fgMask)/(double)fgMask.size().area();

    pMOG2->apply(frame, fgMask,dLearningRate);
    //pKNN->apply(frame, fgMask,dLearningRate);


    //pMOG->apply(frame, fgMaskMOG,dLearningRate);
    //pGMG->apply(frame,fgMaskGMG,dLearningRate);



        //OPENCV 3 MORPHOLOGICAL
    //get the frame number and write it on the current frame
    //erode to get rid to food marks
    //cv::erode(fgMaskMOG2,fgMaskMOG2,kernel, cv::Point(-1,-1),3);
    //Do Close : erode(dilate())
    //cv::morphologyEx(fgMaskMOG2,fgMaskMOG2, cv::MORPH_CLOSE, kernelClose,cv::Point(-1,-1),2);
    //cv::dilate(fgMaskMOG2,fgMaskMOG2,kernel, cv::Point(-1,-1),4);
    //Apply Open Operation dilate(erode())
    //cv::morphologyEx(fgMaskMOG2,fgMaskMOG2, cv::MORPH_OPEN, kernel,cv::Point(-1,-1),2);


//    //Put Info TextOn Frame
//    //Frame Number
//    std::stringstream ss;
//    cv::rectangle(frame, cv::Point(10, 2), cv::Point(100,20),
//              cv::Scalar(255,255,255), -1);
//    ss << nFrame;
//    std::string frameNumberString = ss.str();
//    cv::putText(frame, frameNumberString.c_str(), cv::Point(15, 15),
//            cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));

//    //Count on Original Frame
//    std::stringstream strCount;
//    strCount << "N:" << (nLarva);
//    cv::rectangle(frame, cv::Point(10, 25), cv::Point(100,45), cv::Scalar(255,255,255), -1);
//    cv::putText(frame, strCount.str(), cv::Point(15, 38),
//            cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));

//    char buff[100];
//    //Learning Rate
//    //std::stringstream strLearningRate;
//    std::sprintf(buff,"dL: %0.4f",dLearningRate);
//    //strLearningRate << "dL:" << (double)(dLearningRate);
//    cv::rectangle(frame, cv::Point(10, 50), cv::Point(100,70), cv::Scalar(255,255,255), -1);
//    cv::putText(frame, buff, cv::Point(15, 63),
//            cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));


//    //Time Rate - conv from ms to minutes //TODO: Replace With actual framerate
//    std::sprintf(buff,"t: %0.2fs",nFrame/55.0 );
//    //strTimeElapsed << "" <<  << " m";
//    cv::rectangle(frame, cv::Point(10, 75), cv::Point(100,95), cv::Scalar(255,255,255), -1);
//    cv::putText(frame, buff, cv::Point(15, 88),
//            cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));

//    //Count Fg Pixels // Ratio
//    std::stringstream strFGPxRatio;
//    dblRatioPxChanged = (double)cv::countNonZero(fgMask)/(double)fgMask.size().area();
//    strFGPxRatio << "Dpx:" <<  dblRatioPxChanged;
//    cv::rectangle(frame, cv::Point(10, 100), cv::Point(100,120), cv::Scalar(255,255,255), -1);
//    cv::putText(frame, strFGPxRatio.str(), cv::Point(15, 113),
//            cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));

    return ret; //If False then tell calling function to stop updating
}



//
// Process Larva video, removing BG, detecting moving larva- Setting the learning rate will change the time required
// to remove a pupa from the scene -
//
unsigned int processVideo(cv::Mat& fgMask, MainWindow& window_main, QString videoFilename, QString outFileCSV, unsigned int startFrameCount)
{

    //Speed that stationary objects are removed
    cv::Mat frame,frameMasked,outframe;;
    unsigned int nFrame = startFrameCount; //Current Frame Number

    bPaused =false; //Start Paused

    std::string frameNumberString;

    //?Replicate FG Mask to method specific
    //fgMask.copyTo(fgMaskMOG2);
    //fgMask.copyTo(fgMaskMOG);
    //fgMask.copyTo(fgMaskGMG);


    //Make Variation of FileNames for other Output

    QString trkoutFileCSV = outFileCSV;
    trkoutFileCSV.truncate(trkoutFileCSV.lastIndexOf("."));
    trkoutFileCSV.append("_tracks.csv");
    QString vialcountFileCSV = outFileCSV;
    vialcountFileCSV.truncate(vialcountFileCSV.lastIndexOf("."));
    vialcountFileCSV.append("_N.csv");

    //REPORT
   std::cout << "Tracking data saved to :" << vialcountFileCSV.toStdString()  <<std::endl;
   std::cout << "\t " << outFileCSV.toStdString() <<std::endl;
   std::cout << "\t " << trkoutFileCSV.toStdString()  <<std::endl;




    //create the capture object
    cv::VideoCapture capture(videoFilename.toStdString());
    if(!capture.isOpened())
    {
        //error in opening the video input
        std::cerr << "Unable to open video file: " << videoFilename.toStdString() << std::endl;
        std::exit(EXIT_FAILURE);
    }




    //read input data. ESC or 'q' for quitting
    while( !bExiting && (char)keyboard != 27 )
    {
        try
        {
            //read the current frame
            if(!capture.read(frame))
            {
                if (nFrame == startFrameCount)
                {
                    std::cerr << "Unable to read first frame." << std::endl;
                    nFrame = 0; //Signals To caller that video could not be loaded.
                    exit(EXIT_FAILURE);
                }
                else
                {
                   std::cerr << "Unable to read next frame. So this video Is done." << std::endl;
                   std::cout << nFrame << " frames of Video processed. Move on to next timelapse video? " <<std::endl;
                    ::saveImage(frameNumberString,gstroutDirCSV,frameMasked);
                   //continue;
                   break;
               }
            }
        }catch(const std::exception &e)
        {
            std::cerr << "Error reading frame " << nFrame << "skipping." << std::endl;
            continue;
        }



        //Add frames from Last video
        nFrame = capture.get(CV_CAP_PROP_POS_FRAMES) + startFrameCount;

        //Make Global Roi on 1st frame
        if (nFrame == 1)
        {
            //Add Global Roi - Center - Radius
            ltROI newROI(cv::Point(frame.cols/2,frame.rows/2),ptROI2);
            addROI(newROI);
        }


        processFrame(frame,fgMask,frameMasked,nFrame,outframe);

        window_main.showVideoFrame(outframe,nFrame); //Show On QT Window

        if (bshowMask)
        {
            cv::imshow(gstrwinName + " FG Mask", fgMask);
            //cv::imshow(gstrwinName + " FG Fish Mask", fgMaskFish);
        }


        if (bTracking)
            saveTracks(vfishmodels,trkoutFileCSV,videoFilename,frameNumberString);

        //if (nFrame%10)
       //     keyboard = cv::waitKey( 1 );

        checkPauseRun(&window_main,keyboard,nFrame);


    } //main While loop
    //delete capture object
    capture.release();



    std::cout << "Exiting video processing loop." <<std::endl;

    return nFrame;
}




//Operator for Priority Ordering
bool operator<(const fishModel& a, const fishModel& b)
{
  return a.templateScore > b.templateScore; //Max Heap
}


void UpdateFishModels(cv::Mat& maskedImg_gray,fishModels& vfishmodels,zftblobs& fishblobs,unsigned int nFrame)
{

    qfishModels qfishrank;

    fishModel* pfish = NULL;

    fishModels::iterator ft;

    cv::Size szTempIcon(std::max(fishbodyimg_template.cols,fishbodyimg_template.rows),std::max(fishbodyimg_template.cols,fishbodyimg_template.rows));
    cv::Point rotCentre = cv::Point(szTempIcon.width/2,szTempIcon.height/2);

    cv::Point gptmaxLoc; //point Of Bestr Match

     //Look through Tracks find they have a fish model attached and create if missing
    for (zftblobs::iterator it = fishblobs.begin(); it!=fishblobs.end(); ++it)
    {

        zftblob* fishblob = &(*it);


        ///
        /// Check If Track Centre Point Contains An image that matches a fish template
        ///

        cv::Point centroid = fishblob->pt;
        cv::Point pBound1 = cv::Point(max(0,min(maskedImg_gray.cols,centroid.x-40)), max(0,min(maskedImg_gray.rows,centroid.y-40)));
        cv::Point pBound2 = cv::Point(max(0,min(maskedImg_gray.cols,centroid.x+40)), max(0,min(maskedImg_gray.rows,centroid.y+40)));

        cv::Rect rectFish(pBound1,pBound2);

        cv::rectangle(frameDebugC,rectFish,CV_RGB(20,200,150),2);
        cv::Mat fishRegion(maskedImg_gray,rectFish); //Get Sub Region Image
        double maxMatchScore; //
        int AngleIdx = templatefindFishInImage(fishRegion,gFishTemplateCache,szTempIcon, maxMatchScore, gptmaxLoc,iLastKnownGoodTemplateRow,iLastKnownGoodTemplateCol);


        int bestAngle =AngleIdx*gFishTemplateAngleSteps;
        cv::Point top_left = pBound1+gptmaxLoc;
        cv::Point ptbcentre = top_left + rotCentre;

        bool bModelFound = false;
        //Check Through Models And Find The Closest Fish To This FishBlob
        for ( ft  = vfishmodels.begin(); ft!=vfishmodels.end(); ++ft)
        {
             pfish = ft->second;

             //Check Overlap Of This Model With The Blob - And Whether The Image of this Blob contains something That looks like a fish
             if (pfish->zfishBlob.overlap(pfish->zfishBlob,*fishblob) > 0 && maxMatchScore > gMatchShapeThreshold )
             {
                 //Some existing Fish Can be associated with this Blob - As it Overlaps from previous frame
                bModelFound = true;
                ///Update Model State
                /// But not While Uses Is manually updating/ Modifying Bounding Box (Flags Are set in Mainwindow)
                if (!bStoreThisTemplate && !bDraggingTemplateCentre) //Skip Updating Bound If this round we are saving The Updated Boundary
                    pfish->updateState(fishblob,maxMatchScore,bestAngle,ptbcentre,nFrame);


                //Add To Priority Q So we can Rank
                qfishrank.push(pfish);
             }


        }

       //Check If Fish Was found)
        if (!bModelFound) //Model Does not exist for track - its a new track
        {

            //Make new fish Model
            //fishModel* fish= new fishModel(track,fishblob);
           fishModel* fish= new fishModel(*fishblob);

           fish->updateState(fishblob,maxMatchScore,bestAngle,ptbcentre,nFrame);

           vfishmodels.insert(IDFishModel(fish->ID,fish));
           qfishrank.push(fish);

           std::cout << "New fishmodel: " << fish->ID << " with Template Score :" << fish->templateScore << std::endl;

        }


    }

    ///\brief Make A priority Queue Ranking Candidate Fish with TemplateSCore - Keep Top One Only
    ///     ///Keep Only the Fish with The Max Template Score - Can Add them to priority Queue And just keep top one

    fishModel* pfishBest = 0;
    double maxTemplateScore = 0.0;
    if (qfishrank.size() > 0)
    {
        pfishBest = qfishrank.top();
        maxTemplateScore = pfishBest->templateScore;
    }

    //Delete All FishModels EXCEPT the best Match - Assume 1 Fish In scene / Always Retain 1 Model
    ft = vfishmodels.begin();
    while(ft != vfishmodels.end() && vfishmodels.size() > 1)
    {
        pfish = ft->second;

        if (pfishBest != pfish )
        {
            //Check Ranking Is OK, as long off course that a fishTemplate Has Been Found On This Round -
            //OtherWise Delete The model?
            //Assertion Fails When Old Model Goes Out Of scene and video Is retracked
            //assert(pfish->templateScore < maxTemplateScore || maxTemplateScore == 0);

            std::cout << "Deleted fishmodel: " << pfish->ID << " Low Template Score :" << pfish->templateScore << std::endl;
            ft = vfishmodels.erase(ft);
            delete(pfish);
            continue;
        }
        ++ft; //Increment Iterator
    } //Loop To Delete Other FishModels




}


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
        bTracking = !bTracking;


    if ((char)keyboard == 's')
    {
        std::cout << "Save Image" << endl;
        bSaveImages = !bSaveImages;
        win->saveScreenShot(gstroutDirCSV);

    }

    if ((char)keyboard == 'r')
    {
        std::cout << "Run" << endl;
        bPaused = false;
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

}


void checkPauseRun(MainWindow* win, int keyboard,unsigned int nFrame)
{

//    int ms = 1;
//    struct timespec ts = { ms / 1000, (ms % 1000) * 1000 * 1000 };
//    nanosleep(&ts, NULL);

    QCoreApplication::processEvents(QEventLoop::AllEvents, 1);

        while (bPaused && !bExiting)
        {


            //Wait Until Key to unpause is pressed
            //keyboard = cv::waitKey( 30 );

            QTime dieTime= QTime::currentTime().addSecs(1);
            while (QTime::currentTime() < dieTime)
                QCoreApplication::processEvents(QEventLoop::AllEvents, 100);


            //keyCommandFlag(win,keyboard,nFrame);
        }

}

bool saveImage(std::string frameNumberString,QString dirToSave,cv::Mat& img)
{

    //Save Output BG Masks
    QString imageToSave =   QString::fromStdString( std::string("output_MOG_") + frameNumberString + std::string(".png"));
    //QString dirToSave = qApp->applicationDirPath();

    dirToSave.append("/pics/");
    imageToSave.prepend(dirToSave);

    if (!QDir(dirToSave).exists())
    {
        std::cerr << "Make directory " << dirToSave.toStdString() << std::endl;
        QDir().mkpath(dirToSave);
    }

    bool saved = cv::imwrite(imageToSave.toStdString(), img);
    if(!saved) {
        cv::putText(img,"Failed to Save " + imageToSave.toStdString(), cv::Point(25, 25), cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(250,250,250));
        cv::putText(img,"Failed to Save" + imageToSave.toStdString(), cv::Point(25, 25), cv::FONT_HERSHEY_SIMPLEX, 0.4 , cv::Scalar(0,0,0));
       std::cerr << "Unable to save " << imageToSave.toStdString() << std::endl;
        return false;
    }
    else
    {
     std::cout << "Saved image " << imageToSave.toStdString() <<std::endl;
    }

    //cv::imshow("Saved Frame", img);

    return true;
}


///Don't need this / fish Contours already exist from enhance Fish Mask
int countObjectsviaContours(cv::Mat& srcimg )
{
     cv::Mat imgTraced;
     srcimg.copyTo(imgTraced);
     std::vector< std::vector <cv::Point> > contours; // Vector for storing contour
     std::vector< cv::Vec4i > hierarchy;

     cv::findContours( imgTraced, contours, hierarchy,CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE ); // Find the contours in the image
     for( unsigned int i = 0; i< contours.size(); i=hierarchy[i][0] ) // iterate through each contour.
     {
          cv::Rect r= cv::boundingRect(contours[i]);
          cv::rectangle(imgTraced,r, cv::Scalar(255,0,0),1,8,0);
          cv::rectangle(frameDebugA,r, cv::Scalar(255,0,0),1,8,0);
     }

     //Write text For Count on Original Frame
     std::stringstream strCount;
     strCount << "N:" << ((int)contours.size());

     cv::rectangle(frameDebugA, cv::Point(540, 2), cv::Point(690,20), cv::Scalar(255,255,255), -1);
     cv::putText(frameDebugA, strCount.str(), cv::Point(545, 15),
             cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));

    std::cout << " Larvae  "<< strCount.str() << std::endl;
    //imshow("Contoured Image",frame);


    // To get rid of the smaller object and the outer rectangle created
      //because of the additional mask image we enforce a lower limit on area
      //to remove noise and an upper limit to remove the outer border.

 /* if (contourArea(contours_poly[i])>(mask.rows*mask.cols/10000) && contourArea(contours_poly[i])<mask.rows*mask.cols*0.9){
      boundRect[i] = boundingRect( Mat(contours_poly[i]) );
      minEnclosingCircle( (Mat)contours_poly[i], center[i], radius[i] );
      circle(drawing,center[i], (int)radius[i], Scalar(255,255,255), 2, 8, 0);
      rectangle(drawing,boundRect[i], Scalar(255,255,255),2,8,0);
      num_object++;
}
      */

    return contours.size();
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
    params.minArea = thresh_fishblobarea;
    params.maxArea = 10*thresh_fishblobarea;

    /////An inertia ratio of 0 will yield elongated blobs (closer to lines)
    ///  and an inertia ratio of 1 will yield blobs where the area is more concentrated toward the center (closer to circles).
    params.filterByInertia      = true;
    params.maxInertiaRatio      = 0.8;
    params.minInertiaRatio      = 0.01;


    //params.filterByInertia = true;

    // Set up the detector with default parameters.
    cv::Ptr<cv::SimpleBlobDetector> detector = cv::SimpleBlobDetector::create(params);

    detector->detect( frame, keypoints); //frameMask


    //Mask Is Ignored so Custom Solution Required
    //for (cv::KeyPoint &kp : keypoints)
    ptFishblobs.clear();
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
                     ptFishblobs.push_back(kp);

            //int maskVal=(int)gframeMask.at<uchar>(kp.pt);
            //if (maskVal > 0)
             //keypoints_in_mask.push_back(kp);
        }
    }


    // Draw detected blobs as red circles.
    // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures the size of the circle corresponds to the size of blob
    //frame.copyTo(frameOut,maskimg); //mask Source Image
    cv::drawKeypoints( frameOut, ptFishblobs, frameOut, cv::Scalar(200,20,20), cv::DrawMatchesFlags::DEFAULT );


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
int processFoodBlobs(cv::Mat& frame,cv::Mat& maskimg,cv::Mat& frameOut,std::vector<cv::KeyPoint>& ptFoodblobs)
{

    std::vector<cv::KeyPoint> keypoints;
    //std::vector<cv::KeyPoint> keypoints_in_ROI;
    cv::SimpleBlobDetector::Params params;

    params.filterByCircularity  = false; //a circle has a circularity of 1, circularity of a square is 0.785, and so on.
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
    params.maxArea = 40;

    /////An inertia ratio of 0 will yield elongated blobs (closer to lines)
    ///  and an inertia ratio of 1 will yield blobs where the area is more concentrated toward the center (closer to circles).
    params.filterByInertia      = false;
    params.maxInertiaRatio      = 1.0;
    params.minInertiaRatio      = 0.7;


    //params.filterByInertia = true;

    // Set up the detector with default parameters.
    cv::Ptr<cv::SimpleBlobDetector> detector = cv::SimpleBlobDetector::create(params);

    //\todo - Memory Crash Here
    detector->detect( frame, keypoints); //frameMask


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




int saveTracks(fishModels& vfish,QString filenameCSV,QString filenameVid,std::string frameNumber)
{
    bool bNewFileFlag = true;
    int cnt;
    int Vcnt = 0;

    //Loop Over ROI
    for (ltROIlist::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
    {
        cnt = 1;
        Vcnt++;
        ltROI iroi = (ltROI)(*it);
        //Make ROI dependent File Name
        QFileInfo fiVid(filenameVid);
        QFileInfo fiOut(filenameCSV);
        QString fileVidCoreName = fiVid.completeBaseName();
        QString dirOutPath = fiOut.absolutePath() + "/"; //filenameCSV.left(filenameCSV.lastIndexOf("/")); //Get Output Directory


        char buff[50];
        sprintf(buff,"_tracks_%d.csv",Vcnt);
        dirOutPath.append(fileVidCoreName); //Append Vid Filename To Directory
        dirOutPath.append(buff); //Append extension track and ROI number

        QFile data(dirOutPath);
        if (data.exists())
            bNewFileFlag = false;


        if(data.open(QFile::WriteOnly |QFile::Append))
        {

            QTextStream output(&data);
            if (bNewFileFlag)
                 output << "frameN,fishID,AngleDeg,Centroid_X,Centroid_Y,EyeLDeg,EyeRDeg\n";

            //Save Tracks In ROI
            for (fishModels::iterator it=vfish.begin(); it!=vfish.end(); ++it)
            {
                cnt++;
                fishModel* pfish = it->second;
                cvb::CvLabel cvL = it->first;

                if (iroi.contains(pfish->ptRotCentre))
                    //Printing the position information +
                    //+ lifetime; ///< Indicates how much frames the object has been in scene.
                    //+active; ///< Indicates number of frames that has been active from last inactive period.
                    //+ inactive; ///< Indicates number of frames that has been missing.
                    output << (*pfish) << "\n";
              }
            }
        data.close();

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
            drawROI(frameDebugA);
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
    cv::circle(frameDebugA,ptROI1,3,cv::Scalar(255,0,0),1);
    cv::circle(frameDebugA,ptROI2,3,cv::Scalar(255,0,0),1);

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

void drawROI(cv::Mat& frame)
{
    //frameCpy.copyTo(frame); //Restore Original IMage
    for (std::vector<ltROI>::iterator it = vRoi.begin(); it != vRoi.end(); ++it) {

        ltROI iroi = (ltROI)(*it);
         //cv::rectangle(frame,iroi,cv::Scalar(0,0,250));
         cv::circle(frame,cv::Point(iroi.x() ,iroi.y()),iroi.radius,cv::Scalar(0,0,250),2);

         if (bTracking)
         {
             cv::Point pt1,pt2;
             pt1.x = iroi.centre.x;
             pt1.y = iroi.centre.y;
             pt2.x = pt1.x + iroi.radius;
             pt2.y = pt1.y; //+ iroi.height;

             cv::circle(frame,pt1,3,cv::Scalar(255,0,0),1);
             cv::circle(frame,pt2,3,cv::Scalar(255,0,0),1);


         }
    }
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
/// \brief enhanceFishMask Looks for fish countours and draws them onto the FG mask so as to enhance features
/// This is to recover Background substraction errors -
/// It then uses the fixed image to Find contours *main Internal and External fish contours* using on Masked Image Showing Fish Outline
/// \param frameImg - Raw Input camera input in Mat - colour or gray -
/// \param maskFGImg - Modified Enhanced FG Mask Image
/// \param outFishMask - Mask Enhanced for Fish Blob Detection
/// \param outFoodMaskMask Enhanced for Fish Blob Detection
/// \todo Cross Check Fish Contour With Model Position
/// - Tracker Picks Up Wrong contour Although Template Matching Finds the fish!
///
void enhanceMask(cv::Mat& frameImg, cv::Mat& maskFGImg,cv::Mat& outFishMask,cv::Mat& outFoodMask,std::vector<std::vector<cv::Point> >& outfishbodycontours, std::vector<cv::Vec4i>& outfishbodyhierarchy)
{

int max_thresh = 255;
cv::Mat frameImg_gray;
cv::Mat frameImg_blur;
cv::Mat threshold_output;
//cv::Mat threshold_output_H;
cv::Mat threshold_output_COMB;

//std::vector<std::vector<cv::Point> > fgMaskcontours;
//std::vector<cv::Vec4i> fgMaskhierarchy;



//cv::imshow("MOG2 Mask Raw",maskFGImg);

/////get rid of noise/food marks
////Apply Open Operation dilate(erode())
cv::morphologyEx(maskFGImg,maskFGImg, cv::MORPH_OPEN, kernelOpen,cv::Point(-1,-1),1);
////jOIN bLOB Do Close : erode(dilate())
cv::morphologyEx(maskFGImg,maskFGImg, cv::MORPH_CLOSE, kernelClose,cv::Point(-1,-1),2);



///// Convert image to gray and blur it
cv::cvtColor( frameImg, frameImg_gray, cv::COLOR_BGR2GRAY );
cv::GaussianBlur(frameImg_gray,frameImg_blur,cv::Size(3,3),0);

/// Detect edges using Threshold , A High And  low
g_Segthresh = cv::threshold( frameImg_gray, threshold_output, g_Segthresh, max_thresh, cv::THRESH_BINARY ); // Log Threshold Image + cv::THRESH_OTSU
//cv::adaptiveThreshold(frameImg_gray, threshold_output,max_thresh,cv::ADAPTIVE_THRESH_MEAN_C,cv::THRESH_BINARY,g_Segthresh,0); //Last Param Is const substracted from mean
//ADAPTIVE_THRESH_MEAN_C

//Remove Speckles // Should leave fish INtact
cv::filterSpeckles(maskFGImg,0,3,2 );

/////////////////Make Hollow Mask
//make Inner Fish MAsk /More Accurate Way
//cv::threshold( frameImg_gray, threshold_output_H, g_SegInnerthreshMult * g_Segthresh, max_thresh, cv::THRESH_BINARY ); //Log Threshold Image
//cv::erode(threshold_output,threshold_output_H,kernelOpenfish,cv::Point(-1,-1),g_SegInnerthreshMult);
//Substract Inner from Outer
//cv::bitwise_xor(threshold_output,threshold_output_H,threshold_output_COMB);
///////////////////


//Make Hollow Mask Directly - Broad Approximate -> Grows outer boundary
cv::morphologyEx(threshold_output,threshold_output_COMB, cv::MORPH_GRADIENT, kernelOpenfish,cv::Point(-1,-1),1);

/// Find contours main Internal and External contour using on Masked Image Showing Fish Outline
/// //Used RETR_CCOMP that only considers 1 level children hierachy - I use the 1st child to obtain the body contour of the fish
//First Find What BG Model Considers to be FG
//cv::findContours( maskFGImg, fgMaskcontours,fgMaskhierarchy, cv::RETR_CCOMP,cv::CHAIN_APPROX_TC89_KCOS , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE

outfishbodycontours.clear();
std::vector<std::vector<cv::Point> > fishbodycontours;
std::vector<cv::Vec4i> fishbodyhierarchy;

//Then Use ThresholdImage TO Trace More detailed Contours
cv::findContours( threshold_output_COMB, fishbodycontours,fishbodyhierarchy, cv::RETR_CCOMP,cv::CHAIN_APPROX_TC89_KCOS , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE


outFishMask = cv::Mat::zeros(frameImg_gray.rows,frameImg_gray.cols,CV_8UC1);
threshold_output_COMB.copyTo(outFoodMask);

//std::vector< std::vector<cv::Point> > fishbodyContour_smooth;

///Draw Only the largest contours that should belong to fish
/// \todo Other Match Shapes Could be used here
/// \todo Use WaterShed - Let MOG mask Be FG label and then watershed
int idxFishContour = -1;
std::vector<cv::Point> curve; // THe Fish Contour to use for new Mask
for (int kk=0; kk< fishbodycontours.size();kk++)
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
        int area  = cv::contourArea(fishbodycontours[kk]);

        ///Check Area and then  Find the thresholded Fish Contour std::max(dMeanBlobArea*8,(double)thresh_fishblobarea)
        if (area >  thresh_fishblobarea) //If Contour Is large Enough then Must be fish
        {
            cv::Moments moments =  cv::moments(fishbodycontours[kk]);
            cv::Point centroid;
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
            //std::vector<cv::RotatedRect> rectFeatures;
            //Add Blob To candidate Region of interest Mask
            //idxFishContour = findMatchingContour(fishbodycontours,fishbodyhierarchy,centroid,-1,fgMaskcontours[kk],rectFeatures);
        }

        if (curve.size() >0)
        {
             ///// SMOOTH COntours /////
            double sigma = 1.0;
            int M = round((8.0*sigma+1.0) / 2.0) * 2 - 1; //Gaussian Kernel Size
            assert(M % 2 == 1); //M is an odd number

            //create kernels
            std::vector<double> g,dg,d2g; getGaussianDerivs(sigma,M,g,dg,d2g);

            vector<double> curvex,curvey,smoothx,smoothy,resampledcurveX,resampledcurveY ;
            PolyLineSplit(curve,curvex,curvey);

            std::vector<double> X,XX,Y,YY;
            getdXcurve(curvex,sigma,smoothx,X,XX,g,dg,d2g,false);
            getdXcurve(curvey,sigma,smoothy,Y,YY,g,dg,d2g,false);
            //ResampleCurve(smoothx,smoothy,resampledcurveX,resampledcurveY, 30,false);
            PolyLineMerge(curve,smoothx,smoothy);
            ///////////// END SMOOTHING

            outfishbodycontours.push_back(curve);
            /////COMBINE - DRAW CONTOURS
            //Could Check if fishblob are contained (Doesn't matter if they are updated or not -
            // they should still fall within contour - )
            //cv::drawContours( maskFGImg, fgMaskcontours, kk, CV_RGB(0,0,0), cv::FILLED); //Erase Previous Fish Blob
            cv::drawContours( outFishMask, outfishbodycontours, (int)outfishbodycontours.size()-1, CV_RGB(255,255,255), cv::FILLED); //Draw New One
            cv::drawContours( outFoodMask, outfishbodycontours, (int)outfishbodycontours.size()-1, CV_RGB(0,0,0),3); //Draw New One
      }


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
    }

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
void detectZfishFeatures(cv::Mat& fullImgIn,cv::Mat& fullImgOut, cv::Mat& maskfishFGImg, std::vector<std::vector<cv::Point> >& contours_body,std::vector<cv::Vec4i>& hierarchy_body)
{

    //bool berrorTriangleFit = true;
    //int max_thresh = 255;
    //int idxREyeContour,idxLEyeContour,idxLBodyContour;
    //int idxREyeContourW = -1;
    //int idxLEyeContourW = -1;
    cv::RNG rng(12345);

    cv::Mat maskedImg_gray,maskedfishImg_gray;
    cv::Mat maskfishFeature,maskedfishFeature_blur;


    cv::Mat grad,grad_x, grad_y;
    cv::Mat framelapl,framelapl_buffer, frameCanny;
    // Memory Crash
    std::vector<std::vector<cv::Point> >hull( contours_body.size() );
    std::vector<cv::RotatedRect> rectFeatures; //Fitted Ellipsoids Array

    std::vector<std::vector<cv::Point> > contours_laplace;
    contours_laplace.reserve(contours_body.size());
    std::vector<cv::Vec4i> hierarchy_laplace; //Contour Relationships  [Next, Previous, First_Child, Parent]
    /// Memory Crash on vector Init
    std::vector<std::vector<cv::Point> > contours_laplace_clear; //For contours without markers
    std::vector<cv::Vec4i> hierarchy_laplace_clear; //Contour Relationships  [Next, Previous, First_Child, Parent]
    std::vector<std::vector<cv::Point> > fishfeatureContours( contours_laplace.size() );
    std::vector<cv::RotatedRect> rectfishFeatures; //Fitted Ellipsoids

    std::vector<std::vector<cv::Point> > contours_canny;
    std::vector<cv::Vec4i> hierarchy_canny; //Contour Relationships  [Next, Previous, First_Child, Parent]


    ////// Make Debug Frames ///
    cv::Mat fullImg_colour;
    fullImgIn.convertTo(fullImg_colour,CV_8UC3);
    fullImg_colour.copyTo(frameDebugA);
    fullImg_colour.copyTo(frameDebugB);
    fullImg_colour.copyTo(frameDebugC);
    fullImg_colour.copyTo(frameDebugD);

    //framelapl_buffer.copyTo(framelapl); //Clear Copy On each Iteration


    /// Convert image to gray and blur it
    cv::cvtColor( fullImgIn, maskedImg_gray, cv::COLOR_BGR2GRAY );

    //Make image having masked all fish
    maskedImg_gray.copyTo(maskedfishImg_gray,maskfishFGImg); //Mask The Laplacian //Input Already Masked

    //Blur The Image used to detect  broad features
    cv::GaussianBlur(maskedfishImg_gray,maskedfishFeature_blur,cv::Size(3,3),1,1);

    //cv::Laplacian(maskedfishFeature_blur,framelapl_buffer,CV_8UC1,g_BGthresh);
    //cv::erode(framelapl,framelapl,kernelOpenLaplace,cv::Point(-1,-1),1);
    ///Memory Crash Here Too - remove
    //cv::findContours(framelapl_buffer, contours_laplace_clear,hierarchy_laplace_clear, cv::RETR_CCOMP,cv::CHAIN_APPROX_TC89_L1, cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE
    //cv::imshow("Laplacian Clear",framelapl_buffer);
    cv::Canny( maskedImg_gray, frameCanny, gi_CannyThresSmall,gi_CannyThres  );
    //cv::findContours(frameCanny, contours_canny,hierarchy_canny, cv::RETR_CCOMP,cv::CHAIN_APPROX_NONE , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE


    ////Template Matching Is already Done On Fish Blob/Object
    //Pick The largest dimension and Make A Square
    cv::Size szTempIcon(std::max(fishbodyimg_template.cols,fishbodyimg_template.rows),std::max(fishbodyimg_template.cols,fishbodyimg_template.rows));
   // cv::Point rotCentre = cv::Point(szTempIcon.width/2,szTempIcon.height/2);
    cv::Mat Mrot;

//    ///Detect Head Feature //
//    std::cout << "Match template on #fish:" << vfishmodels.size() << std::endl;
    for (fishModels::iterator it=vfishmodels.begin(); it!=vfishmodels.end(); ++it)
    {
          fishModel* fish = (*it).second;

          //fish->bearingAngle   = AngleIdx;
          if (fish->templateScore < gMatchShapeThreshold)
              continue; //Skip This Model Fish And Check the next one


          //Draw A general Region Where the FIsh Is located, search for template within that region only
          cv::Point centre = fish->ptRotCentre; //top_left + rotCentre;
          //cv::Point centroid = fish->ptRotCentre ; // cv::Point2f(fish->track->centroid.x,fish->track->centroid.y);
          cv::Point pBound1 = cv::Point(max(0,min(maskedImg_gray.cols,centre.x-gFishBoundBoxSize)), max(0,min(maskedImg_gray.rows,centre.y-gFishBoundBoxSize)));
          cv::Point pBound2 = cv::Point(max(0,min(maskedImg_gray.cols,centre.x+gFishBoundBoxSize)), max(0,min(maskedImg_gray.rows,centre.y+gFishBoundBoxSize)));

          cv::Rect rectFish(pBound1,pBound2);

          cv::rectangle(frameDebugC,rectFish,CV_RGB(20,200,150),2); //Identify Fish Region Bound In Cyan Square
          cv::Mat fishRegion(maskedImg_gray,rectFish); //Get Sub Region Image
          //double maxMatchScore; //
          //int AngleIdx = templatefindFishInImage(fishRegion,gFishTemplateCache,szTempIcon, maxMatchScore, gptmaxLoc,iLastKnownGoodTemplateRow,iLastKnownGoodTemplateCol);
           //Check If Fish Was found)
          //fish->templateScore  = maxMatchScore;

          //0 Degrees Is along the Y Axis Looking Upwards
          int bestAngleinDeg = fish->bearingAngle;
          //Set to Global Max Point
         // cv::Point top_left = pBound1+gptmaxLoc;

          ///Write Angle / Show Box

          cv::RotatedRect fishRotAnteriorBox(centre, cv::Size(fishbodyimg_template.cols,fishbodyimg_template.rows),bestAngleinDeg);
          /// Save Anterior Bound
          fish->bodyRotBound = fishRotAnteriorBox;

          stringstream strLbl;
          strLbl << "A: " << bestAngleinDeg;
          cv::putText(fullImgOut,strLbl.str(),fishRotAnteriorBox.boundingRect().br(),CV_FONT_NORMAL,0.4,CV_RGB(250,250,0),1);


          ///Draw a Red Rotated Frame around Detected Body
          cv::Point2f boundBoxPnts[4];
          fishRotAnteriorBox.points(boundBoxPnts);
          for (int j=0; j<4;j++) //Rectangle Body
              cv::line(fullImgOut,boundBoxPnts[j],boundBoxPnts[(j+1)%4] ,CV_RGB(210,00,0),1);

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
          cv::RotatedRect fishEyeBox(ptEyeMid, cv::Size(fishbodyimg_template.cols/2+3,fishbodyimg_template.cols/2+3),bestAngleinDeg);

          // Get Image Region Where the template Match occured
          //- Expand image so as to be able to fit the template When Rotated Orthonormally
          //Custom Bounding Box Needs to allow for RotRect To be rotated Orthonormally
          cv::Rect rectfishAnteriorBound = rectFish; //Use A square // fishRotAnteriorBox.boundingRect();
          cv::Size szFishAnteriorNorm(min(rectfishAnteriorBound.width,rectfishAnteriorBound.height)+4,max(rectfishAnteriorBound.width,rectfishAnteriorBound.height)+4); //Size Of Norm Image
          //Rot Centre Relative To Bounding Box Of UnNormed Image
          cv::Point2f ptFishAnteriorRotCentre = (cv::Point2f)fishRotAnteriorBox.center-(cv::Point2f)rectfishAnteriorBound.tl();

          //Define Regions and Sizes for extracting Orthonormal Fish
          //Top Left Corner of templateSized Rect relative to Rectangle Centered in Normed Img
          cv::Size szTemplateImg = fishbodyimg_template.size();
          //cv::Point ptTopLeftTemplate(szFishAnteriorNorm.width/2-szTemplateImg.width/2,szFishAnteriorNorm.height/2-szTemplateImg.height/2);
          cv::Point ptTopLeftTemplate(rectfishAnteriorBound.width/2-szTemplateImg.width/2,rectfishAnteriorBound.height/2-szTemplateImg.height/2);
          cv::Rect rectFishTemplateBound = cv::Rect(ptTopLeftTemplate,szTemplateImg);
          cv::Size szHeadImg(min(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height),max(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height)*0.75);
//          cv::Point ptTopLeftHead(ptTopLeftTemplate.x,0);//(szFishAnteriorNorm.width/2-szTemplateImg.width/2,szFishAnteriorNorm.height/2-szTemplateImg.height/2);
          cv::Rect rectFishHeadBound = cv::Rect(ptTopLeftTemplate,szHeadImg);


          ///Make Normalized Fish View
           tEllipsoids vell;
           cv::Mat imgTmp, imgFishAnterior,imgFishAnterior_Norm,imgFishHead,imgFishHeadEdge,imgFishHeadProcessed;
           maskedImg_gray.copyTo(imgTmp); //imgTmp Contain full frame Image in Gray
           //Threshold The Match Check Bounds Within Image
           cv::Rect imgBounds(0,0,imgTmp.cols,imgTmp.rows);

           if ( //Looks Like a fish is found, now Check Bounds // gmaxVal > gMatchShapeThreshold &&
               imgBounds.contains(rectfishAnteriorBound.br()) &&
                   imgBounds.contains(rectfishAnteriorBound.tl()))
           {
              imgTmp(rectfishAnteriorBound).copyTo(imgFishAnterior);
              frameCanny(rectfishAnteriorBound).copyTo(imgFishHeadEdge);
              //get Rotated Box Centre Coords relative to the cut-out of the anterior Body - This we use to rotate the image
              ///\note The centre of the Bounding Box could also do




              //cv::Point ptRotCenter = cv::Point(szFishAnteriorNorm.width/2,szFishAnteriorNorm.height/2);
              //cv::Point ptRotCenter = cv::Point(imgFishAnterior.cols/2,imgFishAnterior.rows/2);
              ///Make Rotation MAtrix cv::Point(imgFishAnterior.cols/2,imgFishAnterior.rows/2)
              cv::Point2f ptRotCenter = fishRotAnteriorBox.center - (cv::Point2f)rectfishAnteriorBound.tl();
             // ptRotCenter.x = ptRotCenter.x*cos(bestAngleinDeg*M_PI/180.0);
             // ptRotCenter.y = ptRotCenter.y*sin(bestAngleinDeg*M_PI/180.0);

              cv::Mat Mrot = cv::getRotationMatrix2D( ptRotCenter,bestAngleinDeg,1.0); //Rotate Upwards
              //cv::Mat Mrot = cv::getRotationMatrix2D(-fishRotHeadBox.center,bestAngleinDeg,1.0); //Rotate Upwards

              ///Make Rotation Transformation
              //Need to fix size of Upright/Normed Image
              cv::warpAffine(imgFishAnterior,imgFishAnterior_Norm,Mrot,szFishAnteriorNorm);
              cv::warpAffine(imgFishHeadEdge,imgFishHeadEdge,Mrot,szFishAnteriorNorm);



              ///Store Template Options
              if (bStoreThisTemplate)
              {    //Cut Down To Template Size
                  imgFishAnterior       = imgFishAnterior_Norm(rectFishTemplateBound);
                  addTemplateToCache(imgFishAnterior,gFishTemplateCache,gnumberOfTemplatesInCache);
                  bStoreThisTemplate = false;
              }

              /// Draw Centers for Reference and cleaner Masks
              //Draw  Rotation Centre of Transformation to Norm
              cv::circle(imgFishAnterior,ptRotCenter,3,CV_RGB(100,140,140),1);
              cv::imshow("IsolatedAnterior",imgFishAnterior);

              //Draw Normalized Rotation Centre
              cv::circle(imgFishAnterior_Norm,ptRotCenter,4,CV_RGB(0,0,0),-1);
              imgFishHead           = imgFishAnterior_Norm(rectFishHeadBound);
              //imgFishHead           = imgFishAnterior_Norm;

              //cv::imshow("IsolatedAnteriorTempl",imgFishAnterior);
              //cv::imshow("IsolatedHead",imgFishHead);
              cv::imshow("IsolatedAnteriorNorm",imgFishAnterior_Norm);

              int ret = detectEllipses(imgFishHead,imgFishHeadEdge,imgTmp, bestAngleinDeg,vell,imgFishHeadProcessed);

              if (ret < 2)
                show_histogram("HeadHist",imgFishHead);

              //Paste Eye Processed Head IMage to Into Top Right corner of Larger Image
              cv::Rect rpasteregion(fullImgOut.cols-imgFishHeadProcessed.cols,0,imgFishHeadProcessed.cols,imgFishHeadProcessed.rows );
              imgFishHeadProcessed.copyTo(fullImgOut(rpasteregion));

              ///Print Eye Angle Info
               std::stringstream ss;
               ss.precision(3);
              if (vell.size() > 0)
              {
                  tDetectedEllipsoid lEye = vell.at(0); //L Eye Is pushed 1st
                  fish->leftEye           = lEye;
                  fish->leftEyeTheta      = lEye.rectEllipse.angle;
                  ss << "L:" << fish->leftEyeTheta;
                  cv::putText(fullImgOut,ss.str(),cv::Point(rpasteregion.br().x-75,rpasteregion.br().y+10),CV_FONT_NORMAL,0.4,CV_RGB(250,250,0),1 );
              }

              ss.str(""); //Empty String
              if (vell.size() > 1)
              {
                  tDetectedEllipsoid rEye = vell.at(1); //R Eye Is pushed 2nd
                  fish->rightEye          = rEye;
                  fish->rightEyeTheta     = rEye.rectEllipse.angle;
                  ss << "R:"  << fish->rightEyeTheta;
                  cv::putText(fullImgOut,ss.str(),cv::Point(rpasteregion.br().x-75,rpasteregion.br().y+25),CV_FONT_NORMAL,0.4,CV_RGB(250,250,0),1 );
              }


              ///Do Spine Fitting And Drawing
              if (contours_body.size() > 0)
              {
                fish->fitSpineToContour(maskedImg_gray,contours_body,0,0);
                fish->drawSpine(fullImgOut);
              }

           } //If Fish Img Bound Is With Picture Frame
          ///




    } //For eAch Fish Model

///////////////////
    ///Iterate FISH list - Check If Contour belongs to any fish Otherwise ignore
    for (cvb::CvTracks::const_iterator it = fishtracks.begin(); it!=fishtracks.end(); ++it)
    {
        //Get the  fishmodel associated with this Track
        cvb::CvID trackId = it->first;
        cvb::CvTrack* track = it->second;

     } ///  Track Loop ends here - For Each FishBlob
      ///
        bEyesDetected = false; //Flip Back to off in case it was eye features were marked for saving


//     DEBUG show all contours on Laplace

//        for( size_t i = 0; i< contours_watershed.size(); i++ )
//        {
//             cv::drawContours( frameDebugD, contours_watershed, (int)i, CV_RGB(0,120,250), 1,8,hierarchy_watershed);
//        }

        //cv::drawContours( frameDebugC, contours_watershed, (int)idxLEyeContourW, CV_RGB(220,250,0), 1,8,hierarchy_watershed);
        //cv::drawContours( frameDebugC, contours_watershed, (int)idxREyeContourW, CV_RGB(220,0,0), 1,8,hierarchy_watershed);

        //cv::drawContours( frameDebugC, contours_laplace, idxREyeContour, CV_RGB(60,20,200), 1,8,hierarchy_laplace);
        //cv::drawContours( frameDebugC, contours_laplace, idxLEyeContour, CV_RGB(15,60,210), 1,8,hierarchy_laplace);



    //Draw On Canny Img
    frameCanny.convertTo(frameCanny, CV_8UC3);
    ///DEBUG show all contours -Edge
    for( size_t i = 0; i< contours_canny.size(); i++ )
    {
         //cv::drawContours( frameDebugC, contours_canny, (int)i, CV_RGB(200,0,60), 1,8,hierarchy_canny);
    }



    //cv::imshow("Edges Canny",frameCanny);
    //cv::imshow("Edges Laplace",framelapl);

    //cv::imshow("Debug A",frameDebugA);
    //cv::imshow("Debug B",frameDebugB);
    cv::imshow("Debug C",frameDebugC);
    //cv::imshow("Debug D",frameDebugD);

    //cv::imshow("Threshold H",threshold_output_H);


    //Free Mem
    rectFeatures.clear();
    rectfishFeatures.clear();
}


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



