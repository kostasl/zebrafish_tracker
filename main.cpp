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

#include <QDirIterator>
#include <QDir>
#include <QDebug>


#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
//#include <opencv2/bgsegm.hpp>
#include <opencv2/highgui/highgui.hpp>
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
int g_Segthresh             = 37; //Image Threshold for FIsh Features
int g_SegInnerthreshMult    = 3; //Image Threshold for FIsh Features
int g_BGthresh              = 10; //BG threshold segmentation
int gi_ThresholdMatching    = 10; /// Minimum Score to accept that a contour has been found
bool gOptimizeShapeMatching = false; ///Set to false To disable matchShapes in FindMatching Contour
int gi_CannyThres           = 150;
int gi_CannyThresSmall      = 50; //Aperture size should be odd between 3 and 7 in function Canny
int gi_maxEllipseMajor      = 11; // thres for Hough Transform
int gi_minEllipseMajor      = 7; //thres for Hough Transform
int gi_VotesEllipseThres    = 9; //Votes thres for Hough Transform
int gthresEyeSeg            = 125;
int gnumberOfTemplatesInCache  = 0; //INcreases As new Are Added
const int nTemplatesToLoad = 5; //Number of Templates To Load Into Cache - These need to exist as images in QtResources

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
cvb::CvBlobs blobs; //All Blobs - Updated Ids on everyframe done by cvLabel function
cvb::CvBlobs fishblobs;
cvb::CvBlobs foodblobs;
cvb::CvTracks fishtracks;
cvb::CvTracks foodtracks;
cvb::CvTracks tracks; ///All tracks

//The fish ones are then revaluated using simple thresholding to obtain more accurate contours
fishModels vfishmodels; //Vector containing live fish models

CvFont trackFnt; //Font for Reporting - Tracking

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
    bPaused = true;
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
    //Init Font
    cvInitFont(&trackFnt, CV_FONT_HERSHEY_DUPLEX, 0.4, 0.4, 0, 1);



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
    cvb::cvReleaseBlobs(blobs);
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

    trackVideofiles(window_main);
    //destroy GUI windows
    cv::destroyAllWindows();
    cv::waitKey(0);                                          // Wait for a keystroke in the window

    //pMOG->~BackgroundSubtractor();
    pMOG2->~BackgroundSubtractor();
    //pKNN->~BackgroundSubtractor();
    //pGMG->~BackgroundSubtractor();

    //Empty The Track and blob vectors
    cvb::cvReleaseTracks(tracks);
    cvb::cvReleaseBlobs(blobs);



    std::cout << "Total processing time : mins " << gTimer.elapsed()/60000.0 << std::endl;

    app.quit();

    return app.exec();

}



unsigned int trackVideofiles(MainWindow& window_main)
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

       getBGModelFromVideo(fgMask, window_main,invideoname,outfilename,istartFrame);

       std::cout << "Press r to run Video processing" << std::endl;

       istartFrame = processVideo(fgMask,window_main,invideoname,outfilename,istartFrame);

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
          cv::imshow(gstrwinName + " FG Mask", fgMask);
          //Check For input Control
          keyboard = cv::waitKey( cFrameDelayms );
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
       keyboard = cv::waitKey( cFrameDelayms );
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
           keyboard = cv::waitKey( cFrameDelayms );


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
    pMOG2->apply(outframe, fgMask,dLearningRateNominal);
//    cv::erode(fgMask,fgMask,kernelOpen, cv::Point(-1,-1),1);
    //cv::dilate(fgMaskMOG2,fgMaskMOG2,kernel, cv::Point(-1,-1),4);


    //Draw THe fish Masks more accuratelly by threshold detection - Enhances full fish body detection
    enhanceFishMask(outframe, fgMask,fishbodycontours,fishbodyhierarchy);// Add fish Blobs

    frameMasked = cv::Mat::zeros(frame.rows, frame.cols,CV_8UC3);
    outframe.copyTo(frameMasked,fgMask); //Use Enhanced Mask
    //show the current frame and the fg masks
    cv::imshow(gstrwinName + " FishOnly",frameMasked);


    ///DRAW ROI
    drawROI(outframe);


    lplframe = frameMasked; //Convert to legacy format

    //cvb::CvBlobs blobs;
    ///DO Tracking
    if (bTracking)
    {
       //Simple Solution was to Use Contours To measure Larvae

       // Filters Blobs between fish and food - save into global vectors
        processBlobs(&lplframe,fgMask, blobs,tracks,gstroutDirCSV,frameNumberString,dMeanBlobArea);

        //Here the Track's blob label is updated to the new matching blob
        // Process Food blobs
        cvb::cvUpdateTracks(foodblobs,foodtracks,vRoi, thDistanceFood, inactiveFrameCount,thActive);
        nFood = foodtracks.size();

        // Process Fish blobs
        //ReFilter Let those that belong to fish Contours Detected Earlier
        fishblobs = cvb::cvFilterByContour(fishblobs,fishbodycontours,CV_RGB(10,10,180));
        cvb::cvUpdateTracks(fishblobs,fishtracks,vRoi, thDistanceFish, inactiveFrameCount,thActive);
        nLarva = fishtracks.size();

        //Update Fish Models From Tracks
        UpdateFishModels(vfishmodels,fishtracks);

        //Combine Lists into Tracks before rendering
        tracks.clear();
        //tracks.insert(foodtracks.begin(),foodtracks.end() );
        tracks.insert(fishtracks.begin(),fishtracks.end());

        //
        //saveTracks(tracks,trkoutFileCSV,frameNumberString);

        /// Get Fish Only Image ///

        /// \TODO optimize this. pic is Gray Scale Originally anyway

        detectZfishFeatures(frame,outframe,fgMask,fishbodycontours,fishbodyhierarchy); //Creates & Updates Fish Models

        //Show Tracks
        cvb::cvRenderTracks(tracks, &lplframe, &lplframe,CV_TRACK_RENDER_ID | CV_TRACK_RENDER_PATH,&trackFnt);


        if (bSaveImages)
        {
            saveImage(to_string(nFrame),gstroutDirCSV,frame);
            cv::putText(frameDebugA, "Save ON", cv::Point(15, 600),
                    cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));

        }

    }

    //Save to Disk


    ///

    ///TEXT INFO Put Info TextOn Frame
    //Frame Number
    std::stringstream ss;
    cv::rectangle(outframe, cv::Point(10, 2), cv::Point(100,20),
              cv::Scalar(255,255,255), -1);
    cv::putText(outframe, frameNumberString,  cv::Point(15, 15),
            cv::FONT_HERSHEY_SIMPLEX, 0.4 , cv::Scalar(0,0,0));

    //Count on Original Frame
    std::stringstream strCount;
    strCount << "Nf:" << (nLarva) << " Nr:" << nFood;
    cv::rectangle(outframe, cv::Point(10, 25), cv::Point(120,45), cv::Scalar(255,255,255), -1);
    cv::putText(outframe, strCount.str(), cv::Point(15, 38),
            cv::FONT_HERSHEY_SIMPLEX, 0.4 , cv::Scalar(0,0,0));

    char buff[100];
    //Learning Rate
    //std::stringstream strLearningRate;
    std::sprintf(buff,"dL: %0.4f",dLearningRate);
    //strLearningRate << "dL:" << (double)(dLearningRate);
    cv::rectangle(outframe, cv::Point(10, 50), cv::Point(100,70), cv::Scalar(255,255,255), -1);
    cv::putText(outframe, buff, cv::Point(15, 63),
            cv::FONT_HERSHEY_SIMPLEX, 0.4 , cv::Scalar(0,0,0));

    //Time Rate - conv from ms to minutes

    std::sprintf(buff,"t: %0.2f",gTimer.elapsed()/(1000.0*60.0) );
    //strTimeElapsed << "" <<  << " m";
    cv::rectangle(outframe, cv::Point(10, 75), cv::Point(100,95), cv::Scalar(255,255,255), -1);
    cv::putText(outframe, buff, cv::Point(15, 88),
            cv::FONT_HERSHEY_SIMPLEX, 0.4 , cv::Scalar(0,0,0));

//    //Count Fg Pixels // Ratio
//    std::stringstream strFGPxRatio;
//    dblRatioPxChanged = (double)cv::countNonZero(fgMask)/(double)fgMask.size().area();
//    strFGPxRatio << "Dpx:" <<  dblRatioPxChanged;
//    cv::rectangle(frame, cv::Point(10, 100), cv::Point(100,120), cv::Scalar(255,255,255), -1);
//    cv::putText(frame, strFGPxRatio.str(), cv::Point(15, 113),
//            cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));



}
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

    bPaused =true; //Start Paused

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
            saveTracks(tracks,trkoutFileCSV,frameNumberString);

        keyboard = cv::waitKey( 1 );
        checkPauseRun(&window_main,keyboard,nFrame);


    } //main While loop
    //delete capture object
    capture.release();



    std::cout << "Exiting video processing loop." <<std::endl;

    return nFrame;
}

///
/// \brief UpdateFishModels Create a fish model class attaching a respective fishtrack to it.
///  This prersistence uses the trackid to identify and make informed tracking of fish features across frames
/// \param vfishmodels
/// \param fishtracks
///
/// \note //The whole of  fishModels is deleted after tracking is finished
///
void UpdateFishModels(fishModels& vfishmodels,cvb::CvTracks& fishtracks)
{
    fishModel* pfish = NULL;

     //Look through Tracks find they have a fish model attached and create if missing
    for (cvb::CvTracks::const_iterator it = fishtracks.begin(); it!=fishtracks.end(); ++it)
    {
        cvb::CvTrack* track = it->second;

        fishModels::const_iterator ft =  vfishmodels.find(it->first); //Find model with same Id as the Track - Associated fishModel
        if (track->inactive)
            continue;

        if (ft == vfishmodels.end()) //Model Does not exist for track - its a new track
        {
            //Make Attached FishModel
            cvb::CvBlobs::const_iterator fbt = blobs.find(track->label);
            assert(fbt != blobs.end());
            cvb::CvBlob* fishblob = fbt->second;
            //Make new fish Model
            fishModel* fish= new fishModel(track,fishblob);


            vfishmodels.insert(CvIDFishModel(track->id,fish));

        }
        else ///Some Fish Has that Track ID
        { //Check if pointer is the same Not just the track ID (IDs are re used)
          pfish = ft->second; //Set Pointer to Existing Fish

            //Must point to the same track - OtherWise Replace Fish Model
            if(pfish->track != track)
            {
                delete pfish;
                //Make Attached FishModel
                vfishmodels.erase(track->id); //Replace
                cvb::CvBlobs::const_iterator fbt = blobs.find(track->label);
                assert(fbt != blobs.end());

                cvb::CvBlob* fishblob = fbt->second;
                fishModel* fish= new fishModel(track,fishblob);
                vfishmodels.insert(CvIDFishModel(track->id,fish));


            }

        }



    }

    ///\todo Make A priority Queue Ranking Candidate Fish with TemplateSCore - Keep Top One Only
    double maxTemplateScore = 0; // Save best Templ Score Found Among Fish Models
    //Look Through
    ///Go through Each FishModel And Delete the ones whose tracks are gone
    fishModels::iterator ft = vfishmodels.begin();
    while(ft != vfishmodels.end())
    {
        pfish = ft->second;

        cvb::CvTracks::const_iterator it = fishtracks.find(pfish->ID);

       if (it == fishtracks.end()) //Track No Longer Exists / Delete model
        {
           //ft = vfishmodels.erase(ft); //Only Works On some Compilers
           vfishmodels.erase(ft++);
           std::cout << "Deleted fishmodel: " << pfish->ID << std::endl;
           delete(pfish);
           break;
        }else{ //Track Is inactive Delete Model
           if (pfish->track->inactive)
           {
               std::cout << "Deleted fishmodel: " << pfish->ID << " Track was Inactive t:" << pfish->track->inactive << std::endl;
               //ft = vfishmodels.erase(ft);
               vfishmodels.erase(ft++);
               delete(pfish);
               break;
           }
        }

       //Find Max Template Score Fish
       if (pfish->templateScore > maxTemplateScore)
           maxTemplateScore = pfish->templateScore;

       //Can Only Be reached if above cases eval. false
         ++ft; //Increment Iterator

    }


    ///Keep Only the Fish with The Max Template Score - Can Add them to priority Queue And just keep top one
    ft = vfishmodels.begin();
    while(ft != vfishmodels.end())
    {
        pfish = ft->second;

        if (pfish->templateScore < maxTemplateScore && pfish->templateScore !=0 )
        {
            std::cout << "Deleted fishmodel: " << pfish->ID << " Low Template Score :" << pfish->templateScore << std::endl;
            ft = vfishmodels.erase(ft);
            delete(pfish);
            continue;
        }


        ++ft; //Increment Iterator
    }


} //End Of Update FishModels


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
//             for (fishModels::iterator it=vfishmodels.begin(); it!=vfishmodels.end(); ++it)
//             {
//                 fishModel* fish = (*it).second;
//                   //Let ReleaseTracks Handle This
//                  fish->resetSpine();
//             }
             ReleaseFishModels(vfishmodels);
    }


    //Toggle Show the masked - where blob id really happens
    if ((char)keyboard == 'm')
         bshowMask = !bshowMask;

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
        while (bPaused && !bExiting)
        {
            int ms = 20;
            struct timespec ts = { ms / 1000, (ms % 1000) * 1000 * 1000 };
            nanosleep(&ts, NULL);
            //Wait Until Key to unpause is pressed
            keyboard = cv::waitKey( 30 );

            keyCommandFlag(win,keyboard,nFrame);
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


///
/// \brief processBlobs Separates food from fish using an area filter.
/// It generates a circular mask around the fish so as to allow to process them separatelly
///It renders the food and fish blobs
/// \param srcimg //Legacy format pointer IplImage* to source image
/// \param blobs
/// \param tracks
/// \param outDirCSV
/// \param frameNumberString
/// \param dMeanBlobArea
/// \return
///
int processBlobs(IplImage* srcframeImg,cv::Mat& maskimg,cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString outDirCSV,std::string& frameNumberString,double& dMeanBlobArea)
{

    IplImage  *labelImg;


    ///// Finding the blobs ////////
     int cnt = 0;
     uint minBlobArea = 0;
     uint maxBlobArea = 0;


     IplImage lplfgMaskImg;
///  REGION OF INTEREST - UPDATE - SET
     //*destframeImg        =  srcfullimg; //Convert The Global frame to lplImage
     //cv::Mat fgMaskSurroundFish = cv::Mat::zeros(srcfullimg.rows, srcfullimg.cols,CV_8UC3); //Empty Canvas For Just Fish Blob
     //framefishMaskImg   = (IplImage)fgMaskSurroundFish
     lplfgMaskImg       =  maskimg;

    if (bROIChanged || ptROI2.x != 0)
    {
        //Set fLAG sO FROM now on Region of interest is used and cannot be changed.
        bROIChanged = true;
    }


   //std::cout << "Roi Sz:" << vRoi.size() <<std::endl;
    labelImg=cvCreateImage(cvGetSize(srcframeImg), IPL_DEPTH_LABEL, 1);
    cvb::cvLabel( &lplfgMaskImg, labelImg, blobs );

    cvb::cvFilterByROI(vRoi,blobs); //Remove Blobs Outside ROIs

    cvb::cvBlobAreaStat(blobs,dMeanBlobArea,dVarBlobArea,maxBlobArea,minBlobArea);
    double dsigma = 1.0*std::sqrt(dVarBlobArea);

    ///Separate Fish from Food Blobs
    //copy blobs and then Filter to separate classes
    //Allow only Fish Area Through
    //                                              (CvBlobs &blobs,unsigned int minArea, unsigned int maxArea)
    fishblobs = cvb::cvFilterByArea(blobs,std::min((uint)dMeanBlobArea*8,(uint)thresh_fishblobarea),std::max((uint)(maxBlobArea+dsigma),(uint)thresh_fishblobarea),CV_RGB(10,10,220) ); //Remove Small Blobs

    //Food Blobs filter -> Remove large blobs (Fish)
    ///\todo these blob filters could be elaborated to include moment matching/shape distance
    foodblobs = cvb::cvFilterByArea(blobs,std::max(minBlobArea-dsigma,4.0),(unsigned int)std::min(dMeanBlobArea*3,(double)thresh_fishblobarea/10.0),CV_RGB(0,200,0)); //Remove Large Blobs



    //Debug Show Mean Size Var
    //std::cout << dMeanBlobArea <<  " " << dMeanBlobArea+3*sqrt(dVarBlobArea) <<std::endl;
    ///Go Through Each ROI and Render Blobs - Split Between Fish and Food
    unsigned int RoiID = 0;
    for (std::vector<ltROI>::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
    {
        ltROI iroi = (ltROI)(*it);
        RoiID++;

        //Custom Filtering the blobs for Rendering
        //Count Blobs in ROI
        //Find Fish

        //cvb::CvBlob* fishBlob = cvb::cvLargestBlob(blobs);
        //RENDER FISH
        for (cvb::CvBlobs::const_iterator it = fishblobs.begin(); it!=fishblobs.end(); ++it)
        {
            cvb::CvBlob* blob = it->second;
            cv::Point pnt;
            pnt.x = blob->centroid.x;
            pnt.y = blob->centroid.y;


            if (iroi.contains(pnt))
            {
                //cnt++; //CV_BLOB_RENDER_COLOR
                    cvb::cvRenderBlob(labelImg, blob, &lplfgMaskImg, srcframeImg, CV_BLOB_RENDER_ANGLE | CV_BLOB_RENDER_BOUNDING_BOX, CV_RGB(250,10,10),1);
                    //Make a mask to Surround the fish of an estimated size -  So as to overcome BG Substraction Loses - by redecting countour
                    //cv::circle(fgMaskSurroundFish,cv::Point(blob->centroid.x,blob->centroid.y),((blob->maxx-blob->minx)+(blob->maxy-blob->miny)),CV_RGB(255,255,255),-1);
            }
        }

        //Now Render Food
        for (cvb::CvBlobs::const_iterator it = foodblobs.begin(); it!=foodblobs.end(); ++it)
        {
            cvb::CvBlob* blob = it->second;
            cv::Point pnt;
            pnt.x = blob->centroid.x;
            pnt.y = blob->centroid.y;

            if (iroi.contains(pnt))
            {
                    cvb::cvRenderBlob(labelImg, blob, &lplfgMaskImg, srcframeImg,CV_BLOB_RENDER_COLOR | CV_BLOB_RENDER_CENTROID|CV_BLOB_RENDER_BOUNDING_BOX ,CV_RGB(200,200,0),1.0);
            }
        }

        // render blobs in original image
        //cvb::cvRenderBlobs( labelImg, blobs, &fgMaskImg, &frameImg,CV_BLOB_RENDER_CENTROID|CV_BLOB_RENDER_BOUNDING_BOX | CV_BLOB_RENDER_COLOR);

        //Make File Names For Depending on the Vial - Crude but does the  job
        ///Save Blobs
        if (bSaveBlobsToFile)
        {
            QString strroiFileN = outDirCSV;
            QString strroiFilePos = outDirCSV;
            char buff[150];
            sprintf(buff,"/V%d_pos_N.csv",RoiID);
            strroiFileN.append(buff);
            sprintf(buff,"/V%d_pos.csv",RoiID);
            strroiFilePos.append(buff);

            saveTrackedBlobs(blobs,strroiFilePos,frameNumberString,iroi);
            cnt += saveTrackedBlobsTotals(blobs,tracks,strroiFileN,frameNumberString,iroi);
        }
    } //For Each ROI




    // *always* remember freeing unused IplImages
    cvReleaseImage( &labelImg );

    return cnt;
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


int saveTracks(cvb::CvTracks& tracks,QString filename,std::string frameNumber)
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
        QString strroiFile = filename.left(filename.lastIndexOf("/"));
        char buff[50];
        sprintf(buff,"/V%d_pos_tracks.csv",Vcnt);
        strroiFile.append(buff);

        QFile data(strroiFile);
        if (data.exists())
            bNewFileFlag = false;



        if(data.open(QFile::WriteOnly |QFile::Append))
        {

            QTextStream output(&data);
            if (bNewFileFlag)
                 output << "frameN,TrackID,TrackBlobLabel,Centroid_X,Centroid_Y,Lifetime,Active,Inactive\n";

            //Save Tracks In ROI
            for (cvb::CvTracks::const_iterator it=tracks.begin(); it!=tracks.end(); ++it)
            {
                cnt++;
                cvb::CvTrack* cvT = it->second;
                //cvb::CvLabel cvL = it->first;

                cv::Point pnt;
                pnt.x = cvT->centroid.x;
                pnt.y = cvT->centroid.y;

                if (iroi.contains(pnt))
                    //Printing the position information +
                    //+ lifetime; ///< Indicates how much frames the object has been in scene.
                    //+active; ///< Indicates number of frames that has been active from last inactive period.
                    //+ inactive; ///< Indicates number of frames that has been missing.
                    output << frameNumber.c_str()  << "," << cvT->id  << "," << cvT->label  << "," << cvT->centroid.x << "," << cvT->centroid.y << "," << cvT->lifetime  << "," << cvT->active  << "," << cvT->inactive <<"\n";
              }
            }
        data.close();

   } //Loop ROI
     return cnt;
}
//Mouse Call Back Function
void CallBackFunc(int event, int x, int y, int flags, void* userdata)
{
     if  ( event == cv::EVENT_LBUTTONDOWN )
     {
        bMouseLButtonDown = true;
         //ROI is locked once tracking begins
        if (bPaused && !bROIChanged) //CHANGE ROI Only when Paused and ONCE
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
         //std::cout << "Mouse move over the window - position (" << x << ", " << y << ")" <<std::endl;

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
///
void enhanceFishMask(cv::Mat& frameImg, cv::Mat& maskFGImg,std::vector<std::vector<cv::Point> >& fishbodycontours, std::vector<cv::Vec4i>& fishbodyhierarchy)
{


int max_thresh = 255;
cv::Mat maskfishOnly,frameImg_gray, frameImg_blur,threshold_output,threshold_output_H,threshold_output_COMB;

std::vector<std::vector<cv::Point> > fgMaskcontours;
std::vector<cv::Vec4i> fgMaskhierarchy;

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

//Then Use ThresholdImage TO Trace More detailed Contours
cv::findContours( threshold_output_COMB, fishbodycontours,fishbodyhierarchy, cv::RETR_CCOMP,cv::CHAIN_APPROX_TC89_KCOS , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE


maskfishOnly = cv::Mat::zeros(frameImg_gray.rows,frameImg_gray.cols,CV_8UC1);


std::vector< std::vector<cv::Point> > fishbodyContour_smooth;

///Draw Only the largest contours that should belong to fish
/// \todo Other Match Shapes Could be used here
/// \todo Use WaterShed - Let MOG mask Be FG label and then watershed
int idxFishContour = -1;
std::vector<cv::Point> curve; // THe Fish Contour to use for new Mask
for (int kk=0; kk< fishbodycontours.size();kk++)
{

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
            cv::Point centroid; centroid.x = moments.m10/moments.m00; centroid.y = moments.m01/moments.m00;

            std::vector<cv::RotatedRect> rectFeatures;
            //Add Blob To candidate Region of interest Mask
            //idxFishContour = findMatchingContour(fishbodycontours,fishbodyhierarchy,centroid,-1,fgMaskcontours[kk],rectFeatures);
            curve = fishbodycontours[kk];
        }
        else
        {
            continue; //Next Contour
        }

            ///// SMOOTH COntours /////
            double sigma = 1.0;
            int M = round((3.0*sigma+1.0) / 2.0) * 2 - 1; //Gaussian Kernel Size
            assert(M % 2 == 1); //M is an odd number

            //create kernels
            std::vector<double> g,dg,d2g; getGaussianDerivs(sigma,M,g,dg,d2g);

            vector<double> curvex,curvey,smoothx,smoothy,resampledcurveX,resampledcurveY ;
            PolyLineSplit(curve,curvex,curvey);

            std::vector<double> X,XX,Y,YY;
            getdXcurve(curvex,sigma,smoothx,X,XX,g,dg,d2g,false);
            getdXcurve(curvey,sigma,smoothy,Y,YY,g,dg,d2g,false);
            //ResampleCurve(smoothx,smoothy,resampledcurveX,resampledcurveY, 30,false);
            //PolyLineMerge(curve,smoothx,smoothy);
            PolyLineMerge(curve,smoothx,smoothy);
            fishbodyContour_smooth.push_back(curve);
            ///////////// END SMOOTHING

            /////COMBINE - DRAW CONTOURS
            //Could Check if fishblob are contained (Doesn't matter if they are updated or not -
            // they should still fall within contour - )
            //cv::drawContours( maskFGImg, fgMaskcontours, kk, CV_RGB(0,0,0), cv::FILLED); //Erase Previous Fish Blob
            cv::drawContours( maskfishOnly, fishbodyContour_smooth, (int)fishbodyContour_smooth.size()-1, CV_RGB(255,255,255), cv::FILLED); //Draw New One

            //fishbodycontours[kk].clear();
            //fishbodycontours[kk] = curve;
            //if (idxFishContour > -1)
            fishbodycontours[kk] = curve; //Replace Contour with Smooth Version

}





    //Merge Smoothed Contour Thresholded with BGMAsk //Add the masks so as to enhance fish features
    cv::bitwise_or(maskfishOnly,maskFGImg,maskFGImg);

    //maskfishOnly.copyTo(maskFGImg);

    //threshold_output.copyTo(frameDebugD);

    if (bshowMask)
    {
        cv::imshow("Threshold Out",threshold_output);
        cv::imshow("MOG2 Mask Processed",maskFGImg);
        cv::imshow("Hollow Fish Mask",threshold_output_COMB);
    }
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

    std::vector<std::vector<cv::Point> >hull( contours_body.size() );
    std::vector<cv::RotatedRect> rectFeatures; //Fitted Ellipsoids Array

    std::vector<std::vector<cv::Point> > contours_laplace;
    contours_laplace.reserve(contours_body.size());
    std::vector<cv::Vec4i> hierarchy_laplace; //Contour Relationships  [Next, Previous, First_Child, Parent]
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

    framelapl_buffer.copyTo(framelapl); //Clear Copy On each Iteration



    /// Convert image to gray and blur it
    cv::cvtColor( fullImgIn, maskedImg_gray, cv::COLOR_BGR2GRAY );

    //Make image having masked all fish
    maskedImg_gray.copyTo(maskedfishImg_gray,maskfishFGImg); //Mask The Laplacian //Input Already Masked

    //Blur The Image used to detect  broad features
    cv::GaussianBlur(maskedfishImg_gray,maskedfishFeature_blur,cv::Size(3,3),1,1);

    cv::Laplacian(maskedfishFeature_blur,framelapl_buffer,CV_8UC1,g_BGthresh);
    //cv::erode(framelapl,framelapl,kernelOpenLaplace,cv::Point(-1,-1),1);
    cv::findContours(framelapl_buffer, contours_laplace_clear,hierarchy_laplace_clear, cv::RETR_CCOMP,cv::CHAIN_APPROX_TC89_L1, cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE
    //cv::imshow("Laplacian Clear",framelapl_buffer);
    cv::Canny( maskedImg_gray, frameCanny, gi_CannyThresSmall,gi_CannyThres  );
    //cv::findContours(frameCanny, contours_canny,hierarchy_canny, cv::RETR_CCOMP,cv::CHAIN_APPROX_NONE , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE


////////////USE TEMPLATE MATCHINg /////////////
    cv::Point gptmaxLoc;

    ////No Try Template Matching  Across Angles//
    //Pick The largest dimension and Make A Square
    cv::Size szTempIcon(std::max(fishbodyimg_template.cols,fishbodyimg_template.rows),std::max(fishbodyimg_template.cols,fishbodyimg_template.rows));
    cv::Point rotCentre = cv::Point(szTempIcon.width/2,szTempIcon.height/2);
    cv::Mat Mrot;

//    ///Detect Head Feature //
//    std::cout << "Match template on #fish:" << vfishmodels.size() << std::endl;
    for (fishModels::iterator it=vfishmodels.begin(); it!=vfishmodels.end(); ++it)
    {
          fishModel* fish = (*it).second;
          //Draw A general Region Where the FIsh Is located, search for template within that region only
          cv::Point centroid = cv::Point2f(fish->track->centroid.x,fish->track->centroid.y);
          cv::Point pBound1 = cv::Point(max(0,min(maskedImg_gray.cols,centroid.x-40)), max(0,min(maskedImg_gray.rows,centroid.y-40)));
          cv::Point pBound2 = cv::Point(max(0,min(maskedImg_gray.cols,centroid.x+40)), max(0,min(maskedImg_gray.rows,centroid.y+40)));

          cv::Rect rectFish(pBound1,pBound2);

          cv::rectangle(frameDebugC,rectFish,CV_RGB(20,200,150),2);
          cv::Mat fishRegion(maskedImg_gray,rectFish); //Get Sub Region Image
          double maxMatchScore; //
          int AngleIdx = templatefindFishInImage(fishRegion,gFishTemplateCache,szTempIcon, maxMatchScore, gptmaxLoc,iLastKnownGoodTemplateRow,iLastKnownGoodTemplateCol);
           //Check If Fish Was found)
          fish->templateScore  = maxMatchScore;
          fish->bearingAngle   = AngleIdx;
          if (maxMatchScore < gMatchShapeThreshold)
              continue; //Skip This Model Fish And Check the next one

          //0 Degrees Is along the Y Axis Looking Upwards
          int bestAngleinDeg = AngleIdx*gFishTemplateAngleSteps;
          //Set to Global Max Point
          cv::Point top_left = pBound1+gptmaxLoc;

          ///Write Angle / Show Box
          cv::Point centre = top_left + rotCentre;
          cv::RotatedRect fishRotAnteriorBox(centre, cv::Size(fishbodyimg_template.cols,fishbodyimg_template.rows),bestAngleinDeg);

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
          cv::Rect rectfishAnteriorBound = fishRotAnteriorBox.boundingRect();
          cv::Size szFishAnteriorNorm(min(rectfishAnteriorBound.width,rectfishAnteriorBound.height),max(rectfishAnteriorBound.width,rectfishAnteriorBound.height)); //Size Of Norm Image
          //Rot Centre Relative To Bounding Box Of UnNormed Image
          cv::Point2f ptFishAnteriorRotCentre = (cv::Point2f)fishRotAnteriorBox.center-(cv::Point2f)rectfishAnteriorBound.tl();

          //Define Regions and Sizes for extracting Orthonormal Fish
          //Top Left Corner of templateSized Rect relative to Rectangle Centered in Normed Img
          cv::Size szTemplateImg = fishbodyimg_template.size();
          cv::Point ptTopLeftTemplate(szFishAnteriorNorm.width/2-szTemplateImg.width/2,szFishAnteriorNorm.height/2-szTemplateImg.height/2);
          cv::Rect rectFishTemplateBound = cv::Rect(ptTopLeftTemplate,szTemplateImg);
          cv::Size szHeadImg(min(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height),max(fishRotAnteriorBox.size.width,fishRotAnteriorBox.size.height)*0.75);
          cv::Point ptTopLeftHead(ptTopLeftTemplate.x,0);//(szFishAnteriorNorm.width/2-szTemplateImg.width/2,szFishAnteriorNorm.height/2-szTemplateImg.height/2);
          cv::Rect rectFishHeadBound = cv::Rect(ptTopLeftHead,szHeadImg);


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

              //Make Rotation MAtrix cv::Point(imgFishAnterior.cols/2,imgFishAnterior.rows/2)
              //cv::circle(imgFishAnterior,ptFishAnteriorRotCentre,1,CV_RGB(250,0,0),1);
              cv::imshow("IsolatedAnterior",imgFishAnterior);
              //cv::Point ptRotCenter = cv::Point(szFishAnteriorNorm.width/2,szFishAnteriorNorm.height/2);
              //cv::Point ptRotCenter = cv::Point(imgFishAnterior.cols/2,imgFishAnterior.rows/2);

              cv::Point2f ptRotCenter = fishRotAnteriorBox.center - (cv::Point2f)rectfishAnteriorBound.tl();
              cv::Mat Mrot = cv::getRotationMatrix2D( ptRotCenter,bestAngleinDeg,1.0); //Rotate Upwards
              //cv::Mat Mrot = cv::getRotationMatrix2D(-fishRotHeadBox.center,bestAngleinDeg,1.0); //Rotate Upwards

              ///Make Rotation Transformation
              //Need to fix size of Upright/Normed Image
              cv::warpAffine(imgFishAnterior,imgFishAnterior_Norm,Mrot,szFishAnteriorNorm);
              cv::warpAffine(imgFishHeadEdge,imgFishHeadEdge,Mrot,szFishAnteriorNorm);

              imgFishHead           = imgFishAnterior_Norm(rectFishHeadBound);

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
                  ss << "L:" << lEye.rectEllipse.angle;
                  cv::putText(fullImgOut,ss.str(),cv::Point(rpasteregion.br().x-75,rpasteregion.br().y+10),CV_FONT_NORMAL,0.4,CV_RGB(250,250,0),1 );
              }

              ss.str(""); //Empty String
              if (vell.size() > 1)
              {
                  tDetectedEllipsoid rEye = vell.at(1); //R Eye Is pushed 2nd
                  fish->rightEye          = rEye;
                  ss << "R:"  << rEye.rectEllipse.angle;
                  cv::putText(fullImgOut,ss.str(),cv::Point(rpasteregion.br().x-75,rpasteregion.br().y+25),CV_FONT_NORMAL,0.4,CV_RGB(250,250,0),1 );
              }

              if (bStoreThisTemplate)
              {    //Cut Down To Template Size
                  imgFishAnterior       = imgFishAnterior_Norm(rectFishTemplateBound);
                  addTemplateToCache(imgFishAnterior,gFishTemplateCache,gnumberOfTemplatesInCache);
                  bStoreThisTemplate = false;
              }
           }
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


//void watershedFeatureMethod()
//{
//    //        ////// Draw WATERSHED Labels ////
//    //        /// - \brief With the more accurate positioning of the eye centres now we can obtain
//    //        /// Obtain Accurate FEature Contours for Eyes- Body + Head /////
//    //        ///Make Marked/Labeled Image Using Approx Eye Location
//    //        //This is labelled as uknown Territory with 0
//    //        cv::drawContours( markerEyesImg, contours_body, (int)idxblobContour, CV_RGB(0,0,0), cv::FILLED);
//    //        //markerEyesImg.convertTo(markerEyesImg, CV_32SC1); //CopyTo Changes it To Src Image type?

//    //        //Now Draw Labels on it To Mark L-R Eyes, Head And body region
//    //        cv::Point ptNeck       = pfish->coreTriangle[2]+(pfish->mouthPoint-pfish->coreTriangle[2])*0.3;
//    //        cv::Point ptHead       = pfish->coreTriangle[2]+(pfish->mouthPoint-pfish->coreTriangle[2])*0.5;
//    //        cv::circle(markerEyesImg,pfish->mouthPoint            ,1,CV_RGB(70,70,70),1,cv::FILLED); //Label/Mark Outside Region
//    //        //cv::circle(markerEyesImg,pfish->midEyePoint       ,1,CV_RGB(150,150,150),1,cv::LINE_AA); //Label/Mark Head Region
//    //        cv::circle(markerEyesImg,ptNeck            ,1,CV_RGB(150,150,150),2,cv::LINE_AA); //Label/Mark Head Region
//    //        //cv::circle(markerEyesImg,ptHead       ,1,CV_RGB(150,150,150),1,cv::LINE_AA); //Label/Mark Head Region

//    //        //Can use std::max((int)rectfishFeatures[1].size.width/4,1)
//    //        cv::circle(markerEyesImg,pfish->leftEyePoint,1 ,CV_RGB(255,255,255),1,cv::LINE_AA); //Label/Mark Centre of  Left Eye
//    //        cv::circle(markerEyesImg,pfish->rightEyePoint,1 ,CV_RGB(100,100,100),1,cv::LINE_AA); //Label/Mark Centre Right Eye
//    //        cv::circle(markerEyesImg,pfish->coreTriangle[2],2,CV_RGB(50,50,50),2,cv::LINE_AA); //Label/Mark Body

//    ////        for( size_t i = 0; i< contours_canny.size(); i++ )
//    ////        {
//    ////             cv::drawContours( maskedImg, contours_canny, (int)i, CV_RGB(120,120,120), 1,8,hierarchy_canny);
//    ////        }


//    //        ///END OF  WATERSHED Marking  ///


//            /// WATERSHED Contour Detection ///
//            //markerEyesImg is Input/Ouptu so need to Save Before Marker img is modified in order to debug
//            markerEyesImg.copyTo(tmpMarker,maskfishFeature );
//            //maskedImg.copyTo(imgwatershedShow,maskfishFeature);

//            cv::watershed(fullImg_colour,markerEyesImg); ///Watershed SEGMENTATION

//            //Do Contour on Segmented image

//            //pfish->rightEyeHull.clear();
//            //pfish->leftEyeHull.clear();
//            //cv::findContours(markerEyesImg, contours_watershed,hierarchy_watershed, cv::RETR_CCOMP,cv::CHAIN_APPROX_NONE , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE
//            //idxLEyeContourW = findMatchingContour(contours_watershed,hierarchy_watershed,pfish->leftEyePoint,-1,pfish->leftEyeHull,rectfishFeatures);
//            //idxREyeContourW = findMatchingContour(contours_watershed,hierarchy_watershed,pfish->rightEyePoint,-1,pfish->rightEyeHull,rectfishFeatures);


//            if (idxLEyeContourW!=-1)
//            {
//                cv::convexHull( cv::Mat(contours_watershed[idxLEyeContourW]), pfish->leftEyeHull, false );
//                if (pfish->leftEyeHull.size() > 5)
//                    pfish->leftEyeRect = cv::fitEllipse(pfish->leftEyeHull);
//                else
//                    pfish->leftEyeRect = cv::minAreaRect(pfish->leftEyeHull);
//            }

//            if (idxREyeContourW!=-1)
//            {
//                 cv::convexHull( cv::Mat(contours_watershed[idxREyeContourW]), pfish->rightEyeHull, false );
//                 if (pfish->rightEyeHull.size() > 5)
//                    pfish->rightEyeRect = cv::fitEllipse(pfish->rightEyeHull);
//                 else
//                     pfish->rightEyeRect = cv::minAreaRect(pfish->rightEyeHull);
//            }
//}



/////
///// \brief fitfishCoreTriangle Sets a fixed position to represent fish features Guesses tail point
///// \param maskedfishFeature Image containing a mask of the fish being targeted
///// \param sfish
///// \param contours_body
///// \param idxInnerContour Pass Index for inner fish body contour (as segregated by morph on thresholded image
///// \param idxOuterContour Pass index of the outer whole fish contour
///// \return
/////
//bool fitfishCoreTriangle(cv::Mat& maskfishFeature,cv::Mat& maskedfishImg,fishModel& sfish,std::vector<std::vector<cv::Point> >& contours_body,int idxInnerContour,int idxOuterContour)
//{
//    std::vector<std::vector<cv::Point2f> >triangle; //,triangle_out;
//    bool berrorTriangleFit = false;

//    ///Fit triangle structure to Body
//    // Find Enclosing Triangle of Child contour
//    triangle.resize( contours_body.size() );
//    cv::minEnclosingTriangle(contours_body[idxInnerContour],triangle[idxInnerContour]);
//    cv::minEnclosingTriangle(contours_body[idxOuterContour],triangle[idxOuterContour]);


//    //Check for errors during Fit procedure (they seem to occur on some contours)
//    if (triangle[idxOuterContour].size() > 0 && triangle[idxInnerContour].size() > 0)
//    {
//        //Check All triangle corners
//        for (int k=0;k<3;k++)
//        {
//            //Are coords within bounds?
//            if (triangle[idxOuterContour][k].x <= -10 || triangle[idxInnerContour][k].x <= -10 || triangle[idxOuterContour][k].y <= -10 || triangle[idxInnerContour][k].y <= -10)
//            {
//                berrorTriangleFit = true;
//                break;
//            }
//        }
//    }else
//        berrorTriangleFit = true;


//       if (berrorTriangleFit)
//       {
//           qDebug() << "Error during triangular fit - fitfishCoreTriangle";
//           return berrorTriangleFit; //Exit non critical
//       }




//    //triangle[idxInnerContour] = triangle[idxOuterContour];
////        triangle_out[idxChild] = triangle_out[idxblobContour];


//    ///Map Keypoint Triangle features
//    //Obtain Triangle's side lengths / Find base
//    double dab = cv::norm(triangle[idxOuterContour][0]-triangle[idxOuterContour][1]);
//    double dac = cv::norm(triangle[idxOuterContour][0]-triangle[idxOuterContour][2]);
//    double dbc = cv::norm(triangle[idxOuterContour][1]-triangle[idxOuterContour][2]);

//    cv::Point ptTail;


//    if (dab <= dac && dab <= dbc)
//    {
//        ptTail = triangle[idxOuterContour][2];

//    }
//    else
//    {
//        if (dac <= dab && dac <= dbc)
//        {

//         ptTail  = triangle[idxOuterContour][1];

//        }
//        else
//        { //dbc is the smallest

//          ptTail =   triangle[idxOuterContour][0];
//        }
//    }


//    //Find Inner Triangle Apex - Tail/Body Point
//    dab = cv::norm(ptTail-(cv::Point)triangle[idxInnerContour][0]);
//    dac = cv::norm(ptTail-(cv::Point)triangle[idxInnerContour][1]);
//    dbc = cv::norm(ptTail-(cv::Point)triangle[idxInnerContour][2]);


//    //Find Triangle Width - Set point0 and Point1 to the triangle's base (eyes)
//    if (dab <= dac && dab <= dbc)
//    {
//        sfish.coreTriangle[0] = triangle[idxInnerContour][2];
//        sfish.coreTriangle[1] = triangle[idxInnerContour][1];
//        sfish.coreTriangle[2] = triangle[idxInnerContour][0];

//    }
//    if (dac <= dab && dac <= dbc)
//    {
//        sfish.coreTriangle[0] = triangle[idxInnerContour][0];
//        sfish.coreTriangle[1] = triangle[idxInnerContour][2];
//        sfish.coreTriangle[2] = triangle[idxInnerContour][1];
//    }
//    if (dbc <= dab && dbc <= dac )
//    { //dbc is the smallest
//        sfish.coreTriangle[0] =  triangle[idxInnerContour][0];
//        sfish.coreTriangle[1] =  triangle[idxInnerContour][1];
//        sfish.coreTriangle[2] =  triangle[idxInnerContour][2];
//    }



//    //Set Eye Position
//    //Select Left Right Eye - Set Consistently that point coreTriangle[1] is to the left of [2]
//    cv::Point vecEyeA = sfish.coreTriangle[2] - sfish.coreTriangle[0];
//    cv::Point vecEyeB = sfish.coreTriangle[2] - sfish.coreTriangle[1];

//    //Use As Temp Vars - Gives -Pi  0 +Pi - Convert to 0 2Pi
//    sfish.leftEyeTheta = std::atan2(vecEyeA.y,vecEyeA.x)+M_PI;
//    sfish.rightEyeTheta = std::atan2(vecEyeB.y,vecEyeB.x)+M_PI;

//    if (sfish.leftEyeTheta > sfish.rightEyeTheta )
//    {
//        //use as tmp / Switch R-L eye points Over
//        cv::Point tmp = sfish.coreTriangle[1];
//        sfish.coreTriangle[1] = sfish.coreTriangle[0];
//        sfish.coreTriangle[0] = tmp;

////        double tmpA = sfish.leftEyeTheta ;
////        sfish.leftEyeTheta = sfish.rightEyeTheta;
////        sfish.rightEyeTheta = tmpA;
//    }

//    //Find Position of Body Peak - Relocate Triangle side
//    cv::Point minLoc;
//    cv::Point maxLoc;
//    double minVal,maxVal;
//    //minMaxLoc(InputArray src, double* minVal, double* maxVal=0, Point* minLoc=0, Point* maxLoc=0, InputArray mask=noArray())


//    sfish.tailTopPoint    = ptTail;
//    //if (maxLoc) is contained in triangle?
//    /// \note Problem - Eyes can sometimes be brighter than cyst
//    //cv::minMaxLoc(maskedfishImg,&minVal,&maxVal,&minLoc,&maxLoc,maskfishFeature );
//    //sfish.coreTriangle[2]   = maxLoc; //Now Place index [0] at apex
//    sfish.midEyePoint       = sfish.coreTriangle[0]-(sfish.coreTriangle[0] - sfish.coreTriangle[1])/2;
//    sfish.mouthPoint        = sfish.coreTriangle[2]+(sfish.midEyePoint-sfish.coreTriangle[2])*1.2;

//    ///Temporarly Reposition
//    //sfish.spline[0].x       = sfish.coreTriangle[2].x;
//    //sfish.spline[0].y       = sfish.coreTriangle[2].y;
//    //sfish.calcSpline(sfish.spline);

//    /////DEbug Output
//    ///Draw Fitted inside Triangle
//    if (!berrorTriangleFit)
//    {
//        for (int j=0; j<3;j++)
//           cv::line(frameDebugB,triangle[idxInnerContour][j],triangle[idxInnerContour][(j+1)%3] ,CV_RGB(250,250,00),1,cv::LINE_8);

//        //Draw Fitted outside Triangle
//        for (int j=0; j<3;j++)
//            cv::line(frameDebugB,triangle[idxOuterContour][j],triangle[idxOuterContour][(j+1)%3] ,CV_RGB(255,255,00),1,cv::LINE_8);
//    }

//    //Show Triangle Tail Point on Global Image
//    cv::circle(frameDebugB,ptTail,5,CV_RGB(0,50,200));


//    //Draw body centre point/Tail Top
//    cv::circle(frameDebugB,sfish.coreTriangle[2],5,CV_RGB(20,20,250),1);



//    return berrorTriangleFit; //No error fit
//} //Fit Fish Core Triangle






///// \brief Find point Furthest Along closed outline contour
///// Find point furtherst using shortest paths to each point around a closed contour/outline/spline
//int maxChainDistance(std::vector<cv::Point> vPointChain,int idx,int idy)
//{
//    int antiVertex = idx;
//    int maxminD = 0;
//    //cv::Point ptsrc = vPointChain[idx];
//    int n = vPointChain.size();
//    int dPathL[n]; //Accumulated distance from starting point on Chain Going AntiClockwise
//    int dPathR[n]; //Accumulated distance from starting point on Chain Going Clockwise


//    //Find Antipoint/mirror Points on Chain - where the difference between accumulated distance is minimum
//    dPathL[idx] = 0;
//    dPathR[idx] = 0;

//    int k = idx; //k is index Going In reverse, i going fwd
//    int ringIdx;
//    int ringBIdx;
//    int lastIndexR = k;
//    int lastValL = dPathL[idx];
//    int lastValR = dPathR[idx];

//    for (int i=1;i<n;i++)
//    {

//        k--;

//        ringIdx  = (i+idx)%(n);
//        ringBIdx = (k);

//        //Calculate Leftward And Rightward point Distances - Store in vector
//        dPathL[ringIdx]      = lastValL;
//        dPathR[ringBIdx]     = lastValR;

//        dPathL[ringIdx]   += cv::norm(vPointChain[ringIdx]-vPointChain[(ringIdx+1)%(n)]);
//        dPathR[ringBIdx]  += cv::norm(vPointChain[ringBIdx]-vPointChain[(lastIndexR)%(n)]);

//        lastValL = dPathL[ringIdx];
//        lastValR = dPathR[ringBIdx];
//        lastIndexR = k;
//        if (k==0)
//            k = n; // Do ring Wrap Around
//    }

//    for (int i=0;i<n;i++)
//    {
//        int ringIdx = (i+idx)%(n);
//        //Compare distance to same point from both paths CW CCW and take the shortest one
//        int DistLR = std::min(dPathL[ringIdx],dPathR[ringIdx]);

//        if (DistLR > maxminD)
//        {
//            maxminD = DistLR;
//            //store Index
//            antiVertex = ringIdx;
//        }
//    }

//return antiVertex;
//}



/////
///// \brief findIndexClosesttoPoint - Naive Nearest neighbour finder
///// \param vPointChain array of points
///// \param pt reference point
///// \return index in array of closest point
/////
//int findIndexClosesttoPoint(std::vector<cv::Point> vPointChain,cv::Point pt)
//{
//  int idx;
//  unsigned int minDist  = 0;
//  unsigned int dist     = 0;

//  minDist = 1000000;
//  for (int i=0; i < vPointChain.size();i++)
//  {
//    dist = cv::norm(vPointChain[i]-pt);
//    if (dist < minDist )
//    {
//        idx = i;
//        minDist = dist;
//    }
//  }

//  return idx;
//}

