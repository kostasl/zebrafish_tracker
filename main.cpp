///*
/// \title Zebrafish tracker used in combination with darkfield IR illumination
/// \date Jun 2018
/// \author Konstantinos Lagogiannis
/// \version 1.0
/// \brief Video Analysis software to track zebrafish behaviour from images obtained at high frame rates (>350fps) using darkfield IR
///  illumination(IR light-ring) on a 35mm petridish containing a single animal.
///
/// \note
///     * Chooses input video file, then on the second dialogue choose the text file to export track info in CSV format.
 ///    * The green box defines the region over which the larvae are counted-tracked and recorded to file.
 ///    * Once the video begins to show, use to left mouse clicks to define a new region in the image over which you want to count the larvae.
 ///    * Press p to pause Image. once paused:
 ///    * s to save snapshots in CSV outdir pics subfolder.
 ///    * 2 Left Clicks to define the 2 points of region-of interest for tracking.
 ///    * m to show the masked image of the larva against BG.
 ///    * t Start Tracking
 ///    * f toggle food tracking
 ///    * p to Pause
 ///    * r to UnPause/Run
 ///    * D to delete currently used template from cache
 ///    * R reset fish Spline
 ///    * W Toggle output to CSV file writing
 ///    * T to save current tracked region as new template
 ///    * M Manual Measure Distance (px) -
 ///    * E Manually Set Eye Angles
 ///    * F Manually set prey position (which is then tracked)
 ///    * q Exit Quit application
 ///*
 ///*

///
 ///*  Dependencies : opencv3 (W/O CUDA ) QT5
 ///* /// \details
 ///  Heurestic optimizations:
 ///   * Detection of stopped Larva or loss of features from BG Substraction - via mask correction
 ///   * Filter blobs and maintain separate lists for each class (food/fish)
 ///   * track blobs of different class (food/fish) separatelly so tracks do not interfere
  ///  * Second method of Ellipsoid fitting, using a fast algorithm on edge points
  ///  * Changes template Match region, wide for new blobs, narrow for known fish - Can track at 50fps (06/2018)
  ///  * Combines blob tracking with optic flow at the point of food particle (using Lucas-Kanade) to improve track of prey motion near fish
 ///   * Tail spine is tracking with both, sequential intensity scanning and a variational approach on fitting smoothed fish contour angle and length (estimates fish's tail size)
 ///   * Detect tail and Head points of candidate fish contours: extend tail mask to improve tail spine fitting /Use head pt to inform template matching search region for speed optimizing of larva tracing.

///  \remark OutputFiles
 ///  Data processing:
 ///  * Added Record of Food Count at regular intervals on each video in case, so that even if no fish is being tracked ROI
 ///    the evolution of prey Count in time can be observed. saveTracks outputs a count of prey numbers at a regular interval 1sec, it shows up with fishID 0
 ///
  ///
 /// \bug MOG use under Multi-Processing gives a SegFault in OpenCL - Workaround: Added try block on MOG2, and then flag to switch off OpenCL.
 /// \note Cmd line arguments: /zebraprey~_track --ModelBG=0 --SkipTracked=0  --PolygonROI=1
 ///                           --invideofile=/media/extStore/ExpData/zebrapreyCap/AnalysisSet/AutoSet450fps_18-01-18/AutoSet450fps_18-01-18_WTLiveFed4Roti_3591_009.mp4
 ///                           --outputdir=/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData/TrackerOnHuntEvents_UpTo22Feb/
 ///
 /// \note Example: /zebraprey_track --ModelBG=0 --SkipTracked=0  --PolygonROI=1 --invideolist=VidFilesToProcessSplit1.txt --outputdir=/media/kostasl/Maxtor/KOSTAS/Tracked/
 /// \todo * Add Learning to exclude large detected blobs that fail to be detected as fish - so as to stop fish detection failures
 ///        :added fishdetector class
 ////////

#include <config.h>  // Tracker Constant Defines
#include <larvatrack.h>
#include <ellipse_detect.h>
#include <template_detect.h>
#include <zfttracks.h>
#include <fgmaskprocessing.h>
#include "eyesdetector.h"

#include <errorhandlers.h> // My Custom Mem Fault Handling Functions and Debug

#include <random>

#include <string.h>

#include <cereal/archives/json.hpp> //Data Serialize of EyeDetector
#include <cereal/archives/xml.hpp> //Data Serialize of EyeDetector
#include <fstream>

#include <QDirIterator>
#include <QDir>
#include <QDebug>
//#include <QThread>
#include <QTime>


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

/// CUDA //
/// #include <opencv2/opencv_modules.hpp> //THe Cuda Defines are in here
#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    #include "opencv2/cudaimgproc.hpp"
    #include "opencv2/cudaarithm.hpp"
    #include <opencv2/core/cuda.hpp>
    #include <opencv2/photo/cuda.hpp>
    #include <opencv2/core/cuda_types.hpp>
#endif


#include <GUI/mainwindow.h>
///Curve Smoothing and Matching
#include <CSS/CurveCSS.h>

MainWindow* pwindow_main = nullptr;

// Custom RL optimization of eyeSegmentation and fitting
//init with 20 seg thres states , and 10 eye vergence states
EyesDetector* pRLEye  = new EyesDetector(-5,15,-10,80);    // RL For eye segmentation

//The fish ones are then revaluated using simple thresholding to obtain more accurate contours
fishModels vfishmodels; //Vector containing live fish models
foodModels vfoodmodels;
pointPairs vMeasureLines; //Point pairs defining line distances


int main(int argc, char *argv[])
{
    gTimer.start();

    bROIChanged = true;
    bPaused = false;
    bshowMask = false;
    bTracking = true; //Start By Tracking by default
    bExiting    = false;

    QStringList inVidFileNames; //List of Video Files to Process


    // Get the rdbuf of clog.
    // We will need it to reset the value before exiting.
    auto old_rdbufclog = std::clog.rdbuf();
    auto old_rdbufcerr = std::cerr.rdbuf();


    installErrorHandlers();

    QApplication app(argc, argv);
    //QQmlApplicationEngine engine;

    MainWindow window_main;
    pwindow_main = &window_main;

    window_main.show();

    /// Handle Command Line Parameters //
    const cv::String keys =
        "{help h usage ? |    | print this help  message}"
        "{outputdir   o |    | Dir where To save sequence of images }"
        "{invideofile v |    | Behavioural Video file to analyse }"
        "{invideolist f |    | A text file listing full path to video files to process}"
        "{startframe s | 1  | Video Will start by Skipping to this frame}"
        "{stopframe p | 0  | Video Will stop at this frame}"
        "{startpaused P | 0  | Start tracking Paused On 1st Frame/Need to Run Manually}"
        "{duration d | 0  | Number of frames to Track for starting from start frame}"
        "{logtofile l |    | Filename to save clog stream to }"
        "{ModelBG b | 0  | Learn and Substract Stationary Objects from Foreground mask}"
        "{BGThreshold bgthres | 30  | Absolute grey value used to segment BG (g_Segthresh)}"
        "{SkipTracked t | 0  | Skip Previously Tracked Videos}"
        "{PolygonROI r | 0  | Use pointArray for Custom ROI Region}"
        "{ModelBGOnAllVids a | 1  | Only Update BGModel At start of vid when needed}"
        "{FilterPixelNoise pn | 0  | Filter Pixel Noise During Tracking (Note:This has major perf impact so use only when necessary due to pixel noise. BGProcessing does it by default)}"
        "{DisableOpenCL ocl | 0  | Disabling the use of OPENCL can avoid some SEG faults hit when running multiple trackers in parallel}"
        "{EnableCUDA cuda | 0  | Use CUDA for MOG, and mask processing - if available  }"
        "{HideDataSource srcShow | 0  | Do not reveal datafile source, so user can label data blindly  }"
        "{EyeHistEqualization histEq | 0  | Use hist. equalization to enhance eye detection contrast  }"
        "{TrackFish ft | 1  | Track Fish not just the moving prey }"
        "{MeasureMode M | 0 | Click 2 points to measure distance to prey}"
        ;

    ///Parse Command line Args
    cv::CommandLineParser parser(argc, argv, keys);

    stringstream ssMsg;
    ssMsg<<"Zebrafish Behaviour Tracker"<< std::endl;
    ssMsg<<"--------------------------" << std::endl;
    ssMsg<<"Author : Konstantinos Lagogiannis 2017, King's College London"<<std::endl;
    ssMsg<< "email: costaslag@gmail.com"<<std::endl;
    ssMsg<<"./zebraprey_track <outfolder> <inVideoFile> <startframe=1> <stopframe=0> <duration=inf>"<<std::endl;
    ssMsg<<"(note: output folder is automatically generated when absent)"<<std::endl;
    ssMsg << "Example: \n  Use checkFilesProcessed.sh script to generate list of videos to processes then execute as : " << std::endl;
    ssMsg << "./zebrafish_track -f=VidFilesToProcessSplit1.txt -o=/media/kostasl/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData/Tracked30-11-17/" << std::endl;
    ssMsg << "-Make Sure QT can be found : use export LD_LIBRARY_PATH= path to Qt/5.11.1/gcc_64/lib/  " << std::endl;
    parser.about(ssMsg.str() );

    initGlobalParams(parser,inVidFileNames);

    //If No video Files have been loaded then Give GUI to User //
    if (inVidFileNames.empty())
            inVidFileNames =QFileDialog::getOpenFileNames(nullptr, "Select videos to Process",gstrinDirVid.toStdString().c_str(),
                                                          "Video file (*.mpg *.avi *.mp4 *.h264 *.mkv *.tiff *.png *.jpg *.pgm)", nullptr, nullptr);

    // get the applications dir path and expose it to QML
    //engine.load(QUrl(QStringLiteral("qrc:///main.qml")));


#ifdef    _ZTFDEBUG_
    cv::namedWindow("Debug D",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    cv::namedWindow("Debug A",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
    cv::namedWindow("Debug B",CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);


    frameDebugA = cv::Mat::zeros(640, 480, CV_8U);
    frameDebugB = cv::Mat::zeros(640, 480, CV_8U);
    frameDebugD = cv::Mat::zeros(640, 480, CV_8U);
#endif

    frameDebugC = cv::Mat::zeros(640, 480, CV_8U);
    initROI();

    /// create Background Subtractor objects
    //(int history=500, double varThreshold=16, bool detectShadows=true

    //Init MOG BG substractor
    initBGSubstraction();

    int iLoadedTemplates = initDetectionTemplates();
    pwindow_main->nFrame = 1;
    pwindow_main->LogEvent(QString::number(iLoadedTemplates) + QString("# Templates Loaded "));



   /// Start Tracking of Video Files ///
   try{
        //app.exec();
        std::clog << gTimer.elapsed()/60000.0 << " >>> Start frame: " << uiStartFrame << " StopFrame: " << uiStopFrame << " <<<<<<<<<"  << std::endl;
        trackVideofiles(window_main,gstroutDirCSV,inVidFileNames,uiStartFrame,uiStopFrame);

    }catch (char *e)
    {
        //printf("Exception Caught: %s\n",e);
        qDebug() << "[Error] >>> Exception Caught while processing: " << outfishdatafile.fileName();
        std::cerr << "[Error] Memory Allocation Error :" << e;
        //std::cerr << "Memory Allocation Error! - Exiting";
        std::cerr << "[Error] Close And Delete Current output file: " << outfishdatafile.fileName().toStdString() ;
        closeDataFile(outfishdatafile);
        removeDataFile(outfishdatafile);
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
    kernelDilateMOGMask.release();
    kernelOpen.release();
    gLastfishimg_template.release();

    gFishTemplateCache.release();

    // Save State Space of Reinforcement Learning
    pRLEye->SaveState();
    window_main.LogEvent("[INFO] Saved EyeDetector State.");
    delete pRLEye;//Destroy EyeSeg Assistant



    //gFishTemplateCache.deallocate();

    //app.quit();
    window_main.close();
    cv::destroyAllWindows();


    // Reset the rdbuf of clog.
     std::clog.rdbuf(old_rdbufclog);
     std::cerr.rdbuf(old_rdbufcerr);


    app.quit();
    //Catch Any Mem Alloc Error
    ///\note ever since I converted gFishCache to UMat, a deallocation error Is Hit - UMat was then Removed
    /// This Is  KNown But When OpenCL Is False https://github.com/opencv/opencv/issues/8693
    std::exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;

}



unsigned int trackVideofiles(MainWindow& window_main,QString outputFileName,QStringList invideonames,unsigned int istartFrame = 0,unsigned int istopFrame = 0)
{
    cv::Mat fgMask;
    cv::Mat bgStaticMask;

    QString invideoname = "*.mp4";
    QString nextvideoname;
    //Show Video list to process
    //std::cout << "Video List To process:" <<std::endl;
    if (!bBlindSourceTracking)
    {
        window_main.LogEvent("Video List To process:");
        for (int i = 0; i<invideonames.size(); ++i)
        {
           invideoname = invideonames.at(i);
           //std::cout << "*" <<  invideoname.toStdString() << std::endl;
           window_main.LogEvent(invideoname );
        }
    }

    //Go through Each Image/Video - Hold Last Frame N , make it the start of the next vid.
    for (int i = 0; i<invideonames.size() && !bExiting; ++i)
    {

       //Empty Vector of Fish Models - and Reset ID Counter // gi_MaxFoodID = gi_MaxFishID = 1; - Done In Release
       ReleaseFishModels(vfishmodels);
       ReleaseFoodModels(vfoodmodels);

       invideoname = invideonames.at(i);

       nextvideoname = invideonames.at(std::min(invideonames.size()-1,i+1));
       gstrvidFilename = invideoname; //Global

       std::clog << gTimer.elapsed()/60000.0 << " Now Processing : "<< invideoname.toStdString() << " StartFrame: " << istartFrame << std::endl;
       //cv::displayOverlay(gstrwinName,"file:" + invideoname.toStdString(), 10000 );

       ///Open Output File Check If We Skip Processed Files
       if ( !openDataFile(outputFileName,invideoname,outfishdatafile) )
       {
            if (bSkipExisting) //Failed Due to Skip Flag
                 continue; //Do Next File
       }else
           writeFishDataCSVHeader(outfishdatafile);

       ///Open Output File Check If We Skip Processed Files
       if (openDataFile(outputFileName,invideoname,outfooddatafile,"_food") )
           writeFoodDataCSVHeader(outfooddatafile);
       else
           pwindow_main->LogEvent("[Error] Cannot open tracked prey data file.");


       // Removed If MOG Is not being Used Currently - Remember to Enable usage in enhanceMask if needed//
       if ((bUseBGModelling && gbUpdateBGModel) || (bUseBGModelling && gbUpdateBGModelOnAllVids) )
       {
           //If BG Model Returns >1 frames
            if (getBGModelFromVideo(bgStaticMask, window_main,invideoname,outfilename,MOGhistory))
            {
                cv::dilate(bgStaticMask,bgStaticMask,kernelDilateMOGMask,cv::Point(-1,-1),2);
                cv::morphologyEx(bgStaticMask,bgStaticMask,cv::MORPH_CLOSE,kernelDilateMOGMask,cv::Point(-1,-1),4); //
                cv::bitwise_not ( bgStaticMask, bgStaticMask ); //Invert Accumulated MAsk TO Make it an Fg Mask


                //Next Video File Most Likely belongs to the same Experiment / So Do not Recalc the BG Model
                if (compString(invideoname,nextvideoname) < 3 && !gbUpdateBGModelOnAllVids)
                    gbUpdateBGModel = false; //Turn Off BG Updates
            }

       }

        //Next File Is Different Experiment, Update The BG
       if (compString(invideoname,nextvideoname) > 2  )
           gbUpdateBGModel = true;

       QFileInfo fiVidFile(invideoname);
       if (bBlindSourceTracking)
          window_main.setWindowTitle("Labelling Hunt Event");
       else
         window_main.setWindowTitle("Tracking:" + fiVidFile.completeBaseName() );
       window_main.nFrame = 0;
       window_main.tickProgress(); //Update Slider



       //if (bStaticAccumulatedBGMaskRemove) //Hide The PopUp
        //cv::destroyWindow("Accumulated BG Model");

        //Can Return 0 If DataFile Already Exists and bSkipExisting is true
        uint ret = processVideo(bgStaticMask,window_main,invideoname,outfishdatafile,istartFrame,istopFrame);

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


/// \brief
///
void processFrame(MainWindow& window_main,const cv::Mat& frame,cv::Mat& bgStaticMask, unsigned int nFrame,cv::Mat& outframe,cv::Mat& outframeHeadEyeDetected,cv::Mat& frameHead)
{
    cv::Mat frame_gray,fgMask,fgFishMask,fgFishImgMasked;
    cv::Mat fgFoodMask,bgROIMask;


    //std::vector<cv::KeyPoint>  ptFoodblobs;
    //std::vector<cv::KeyPoint> ptFishblobs;
    zftblobs ptFoodblobs;
    zftblobs ptFishblobs;

    std::vector<std::vector<cv::Point> > fishbodycontours;
    std::vector<cv::Vec4i> fishbodyhierarchy;

    unsigned int nLarva         =  0;
    unsigned int nFood          =  0;
    double dblRatioPxChanged    =  0.0;

    QString frameNumberString;
    frameNumberString = QString::number(nFrame);

    assert(!frame.empty());

    //For Morphological Filter
    ////cv::Size sz = cv::Size(3,3);
    //frame.copyTo(inputframe); //Keep Original Before Painting anything on it
    //update the background model
    //OPEN CV 2.4
    // dLearningRate is now Nominal value

     // frame.copyTo(outframe); //Make Replicate On which we draw output


    ///DRAW ROI
    if (bRenderToDisplay)
        drawAllROI(outframe);

    if (bMeasure2pDistance)
        drawUserDefinedPoints(outframe);
    //Redraw ROI Mask
    if (bROIChanged)
    {
        bgROIMask = cv::Mat::zeros(frame.rows,frame.cols,CV_8UC1);
        vRoi.at(0).drawMask(bgROIMask);
    }


    //lplframe = frameMasked; //Convert to legacy format

    //cvb::CvBlobs blobs;
    ///DO Tracking
    if (bTracking)
    {
       //Simple Solution was to Use Contours To measure LUarvae
        //cvtColo frame_grey
        //Draw THe fish Masks more accuratelly by threshold detection - Enhances full fish body detection
    //    enhanceFishMask(outframe, fgMask,fishbodycontours,fishbodyhierarchy);// Add fish Blobs
        if (frame.channels() > 2)
            cv::cvtColor( frame, frame_gray, cv::COLOR_BGR2GRAY);
        else
            frame.copyTo(frame_gray);
        // Save COpy as Last Frame
        gframeCurrent.copyTo(gframeLast);
        frame_gray.copyTo(gframeCurrent); //Copy To global Frame



        /// DO BG-FG SEGMENTATION MASKING and processing///
        /// \brief processMasks
#if _DEBUG
        cv::imshow("FG Image", frame_gray - gframeBGImage);
#endif
        processMasks(frame_gray,bgStaticMask,fgMask,dLearningRateNominal); //Applies MOG if bUseBGModelling is on
        enhanceMask(frame_gray,fgMask,fgFishMask,fgFoodMask,fishbodycontours, fishbodyhierarchy); //Generates separate masks for Fish/Prey and Draws Fish Contourmask

        //Combine Roi Mask Only For The foodMask
        cv::bitwise_and(bgROIMask,fgFoodMask,fgFoodMask);


        /// Choose FG image prior to template matching
        if (bApplyFishMaskBeforeFeatureDetection)
            frame_gray.copyTo(fgFishImgMasked,fgMask); //fgMask / NOT fish Mask //Use Enhanced Mask
        else
            frame_gray.copyTo(fgFishImgMasked); //Full Frame Image Used


        cv::Mat maskedImg_gray;
        /// Convert image to gray and blur it
        cv::cvtColor( frame, maskedImg_gray, cv::COLOR_BGR2GRAY );


        //Can Use Fish Masked fgFishImgMasked - But Templates Dont Include The masking
        processFishBlobs(fgFishImgMasked,fgFishMask, outframe , ptFishblobs);

        ///Update Fish Models Against Image and Tracks - Obtain Bearing Angle Using Template
        //Can Use Fish Masked - But Templates Dont Include The masking
        //UpdateFishModels(fgFishImgMasked,vfishmodels,ptFishblobs,nFrame,outframe);
        if (bTrackFish)
        {
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
        }

        nLarva = vfishmodels.size();

        ///\todo Keep A Global List of all tracks?

        /// Isolate Head, Get Eye models, and Get and draw Spine model
        if (nLarva > 0 &&  (bTrackFish))
            //An Image Of the Full Fish Is best In this Case
            //Do Not Use Masked Fish Image For Spine Fitting
            detectZfishFeatures(window_main, frame_gray,outframe,frameHead,outframeHeadEyeDetected,fgFishImgMasked,fishbodycontours,fishbodyhierarchy); //Creates & Updates Fish Models

        ///////  Process Food Blobs ////
        // Process Food blobs

        if (bTrackFood)
        {
            processFoodBlobs(frame_gray,fgFoodMask, outframe , ptFoodblobs); //Use Just The Mask
            UpdateFoodModels(maskedImg_gray,vfoodmodels,ptFoodblobs,nFrame,true); //Make New Food Models based on identified Blob

            if (nFrame > gcMinFoodModelActiveFrames)
            {
                processFoodOpticFlow(frame_gray, gframeLast ,vfoodmodels,nFrame,ptFoodblobs ); // Use Optic Flow
                UpdateFoodModels(maskedImg_gray,vfoodmodels,ptFoodblobs,nFrame,false); //Update but no new Food models
            }

            //else
            //cv::imshow("Food Mask",fgFoodMask); //Hollow Blobs For Detecting Food

            //cv::drawKeypoints(outframe,ptFoodblobs)
            if (ptFoodblobs.size() >0)
                cv::drawKeypoints( outframe, ptFoodblobs, outframe, cv::Scalar(20,70,255,60), cv::DrawMatchesFlags::DRAW_RICH_KEYPOINTS );




            ///Draw Food Tracks
            foodModels::iterator ft = vfoodmodels.begin();
            nFood = 0;
            while (ft != vfoodmodels.end() && bRenderToDisplay)
            {

                foodModel* pfood = ft->second;
                assert(pfood);

                // Render Food that has been on for A Min of Active frames / Skip unstable Detected Food Blob - Except If Food is being Tracked
                if ( (pfood->isNew) && (!pfood->isTargeted))
                {
                    ++ft; //Item Is not Counted
                    continue;
                }

                if (pfood->isTargeted) //Draw Track Only on Targetted Food
                    zftRenderTrack(pfood->zTrack, frame, outframe,CV_TRACK_RENDER_ID | CV_TRACK_RENDER_HIGHLIGHT  | CV_TRACK_RENDER_PATH | CV_TRACK_RENDER_BOUNDING_CIRCLE, CV_FONT_HERSHEY_PLAIN, trackFntScale*1.2 ); //| CV_TRACK_RENDER_BOUNDING_BOX

                else{
                if (pfood->isActive)
                {
                        nFood++; //only count the rendered Food Items ie. Active Ones
                        zftRenderTrack(pfood->zTrack, frame, outframe,CV_TRACK_RENDER_ID | CV_TRACK_RENDER_BOUNDING_CIRCLE , CV_FONT_HERSHEY_PLAIN,trackFntScale );
                }
                    else
                        zftRenderTrack(pfood->zTrack, frame, outframe,CV_TRACK_RENDER_ID | CV_TRACK_RENDER_BOUNDING_BOX , CV_FONT_HERSHEY_PLAIN,trackFntScale );
                 }
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
    {
        drawFrameText(window_main,nFrame,nLarva,nFood,outframe);

    }

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
               CV_RGB(10,10,10), CV_FILLED,LINE_8);
    cv::putText(outframe, frameNumberString.toStdString(),  cv::Point(15, 15),
            trackFnt, trackFntScale ,  CV_RGB(150,80,50));

    //Count on Original Frame
    std::stringstream strCount;
    strCount << "Nf:" << (nLarva) << " Nr:" << nFood;
    cv::rectangle(outframe, cv::Point(10, 25), cv::Point(80,45),  CV_RGB(10,10,10), CV_FILLED);
    cv::putText(outframe, strCount.str(), cv::Point(15, 38),
           trackFnt, trackFntScale ,  CV_RGB(150,80,50));

/*
 *     //Report Time
    std::sprintf(buff,"t: %0.2f",gTimer.elapsed()/(1000.0*60.0) );
    //strLearningRate << "dL:" << (double)(dLearningRate);
    cv::rectangle(outframe, cv::Point(10, 50), cv::Point(50,70), cv::Scalar(10,10,10), -1);
    cv::putText(outframe, buff, cv::Point(15, 63),
            trackFnt, trackFntScale , CV_RGB(250,250,0));
*/

} //DrawFrameText



//
// Process Larva video, removing BG, detecting moving larva- Setting the learning rate will change the time required
// to remove a pupa from the scene -
//
unsigned int processVideo(cv::Mat& bgStaticMask, MainWindow& window_main, QString videoFilename, QFile& outdatafile, unsigned int startFrameCount,unsigned int stopFrame=0)
{

    QElapsedTimer otLastUpdate; //Time Since Last Progress Report
    otLastUpdate.start();
    //Speed that stationary objects are removed
    cv::Mat frame,outframe,outframeHeadEyeDetect,outframeHead; //bgROIMask,bgMaskWithRoi
    outframeHead                = cv::Mat::zeros(gszTemplateImg.height,gszTemplateImg.width,CV_8UC1); //Initiatialize to Avoid SegFaults

    unsigned int nFrame         = 0;
    unsigned int nErrorFrames   = 0;

    QString frameNumberString;
    //?Replicate FG Mask to method specific
    //fgMask.copyTo(fgMaskMOG2);
    //fgMask.copyTo(fgMaskMOG);
    //fgMask.copyTo(fgMaskGMG);

    bPaused = false;
    //Make Variation of FileNames for other Output

    //QString trkoutFileCSV = outFileCSV;

    //Make Global Roi on 1st frame if it doesn't prexist
    //Check If FG Mask Has Been Created - And Make A new One



    //create the capture object
    cv::VideoCapture capture(videoFilename.toStdString());
    if(!capture.isOpened())
    {
        //error in opening the video input
        window_main.LogEvent("[ERROR] Failed to open video capture device");
        std::cerr << gTimer.elapsed()/60000.0 << " [Error] Unable to open video file: " << videoFilename.toStdString() << std::endl;
        return 0;
        //std::exit(EXIT_FAILURE);
    }


    gfVidfps  = capture.get(CAP_PROP_FPS);

    gcMaxFoodModelInactiveFrames  = gfVidfps*2; //Number of frames inactive (Not Matched to a Blob) until track is deleted
    gcMinFoodModelActiveFrames    = gfVidfps/5;
    gFoodReportInterval           = gfVidfps; //Report Food every second

    uint totFrames = capture.get(CV_CAP_PROP_FRAME_COUNT);
    window_main.setTotalFrames(totFrames);
    //window_main.nFrame = nFrame;

    //  Check If it contains no Frames And Exit
    if (totFrames < 2)
    {
        window_main.LogEvent("[ERROR] This Video File is empty ");
        capture.release();
        return 0;
    }

    if (!bBlindSourceTracking)
    {
        QFileInfo vidFileInfo(videoFilename);
        window_main.LogEvent(" **Begin Processing: " + vidFileInfo.completeBaseName());
        std::cout << " **Begin Processing: " << vidFileInfo.completeBaseName().toStdString() << std::endl; //Show Vid Name To StdOUt
     }else
        window_main.LogEvent("** Begin Processing of video file ** ");

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
    frameNumberString = QString("%1").arg(nFrame, 5, 10, QChar('0')); //QString::number(nFrame);
    window_main.nFrame = nFrame;

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

        if (!bPaused  )
        {
            nFrame = capture.get(CV_CAP_PROP_POS_FRAMES);
            window_main.nFrame = nFrame; //Update The Frame Value Stored in Tracker Window
            window_main.tickProgress();
        }


        if (nFrame == startFrameCount && !bPaused) //Only Switch Tracking On When Running Vid.
        {
            bTracking = true;
        }

         frameNumberString = QString("%1").arg(nFrame, 5, 10, QChar('0')); //QString::number(nFrame); //QString::number(nFrame); //Update Display String Holding FrameNumber

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

    //Pause on 1st Frame If Flag Start Paused is set
    if (bStartPaused && nFrame == startFrameCount && !bPaused)
    {
        bPaused =true; //Start Paused //Now Controlled By bstartPaused
        pwindow_main->LogEvent(QString("[info]>> Video Paused<<"));
    }

    nErrorFrames = 0;

    //If No Mask Exist Then Make A blank Full On One
    if (bgStaticMask.cols == 0)  {
       bgStaticMask = cv::Mat::ones(frame.rows,frame.cols,CV_8UC1);
    }

    //Blank Drawing Canvas for output - We then Blend with original Image
    if (bRenderWithAlpha)
        outframe = cv::Mat::zeros(frame.rows,frame.cols,frame.type());
    else
        outframe = frame.clone();


        //Pass Processed bgMask which Is then passed on to enhanceMask
    processFrame(window_main,frame,bgStaticMask,nFrame,outframe,outframeHeadEyeDetect,outframeHead);

        double alpha = 0.5;
        if (bRenderToDisplay)
        {
            //Simulated Alpha Channels Causes delays!
            if (bRenderWithAlpha)
                cv::addWeighted(frame,1.0,outframe,1.0-alpha,0.0,outframe);

            ///Paste Eye Processed Head IMage to Into Top Right corner of Larger Image
            rect_pasteregion.width = outframeHeadEyeDetect.cols;
            rect_pasteregion.height = outframeHeadEyeDetect.rows;
            if (outframeHeadEyeDetect.u)
              outframeHeadEyeDetect.copyTo(outframe(rect_pasteregion) ) ;

//          cv::imshow("headDetect",outframeHeadEyeDetect);

            window_main.showVideoFrame(outframe,nFrame); //Show On QT Window
            window_main.showInsetimg(outframeHead); //Shows the edge processing img used for eye detection
        }

        // Switch Off Render To Display If In Offline Tracking Mode
        //-Placed Here to Allow 1Frame Short Toggles Of Rendering
        bRenderToDisplay = !bOffLineTracking;
        //frame.copyTo(frameDebugD);
        //cv::imshow("Debug D",frameDebugD);

        /// Report Memory Usage Periodically - Every realtime Second//
        if (!bPaused && (nFrame % (uint)gfVidfps) == 0 || !bPaused && nFrame == 2)
        {
            double rss,vm;
            process_mem_usage(vm, rss);
            //std::clog << "Delta Memory VM: " << vm/1024.0 << "MB; RSS: " << rss/1024.0 << "MB" << std::endl;
            //Show Memory Consumption

            std::stringstream ss;
            ss.precision(4);
            float fFps = gfVidfps/((float)otLastUpdate.elapsed()/1000.0);
            ss  << " [Progress] " << nFrame <<"/" << totFrames << " fps:" <<  fFps << " D Memory VM: " << vm/1024.0 << "MB; RSS: " << rss/1024.0 << "MB";
            window_main.LogEvent(QString::fromStdString(ss.str()));

            otLastUpdate.restart();
            // Render Next Frame To Display
            if (bOffLineTracking)
                bRenderToDisplay = true;

            ss.str(std::string()); //Clear
          //Report MOG Mixtures
            ss << "[Progress] MOGMixtures : " << pMOG2->getNMixtures() << " VarThres:" << pMOG2->getNMixtures() << " VarGen:" << pMOG2->getVarThresholdGen();
            window_main.LogEvent(QString::fromStdString(ss.str()));

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
        if (bTracking && !bPaused && bRecordToFile)
        {
            if (!saveTracks(vfishmodels,vfoodmodels,outfishdatafile,frameNumberString))
            {
                bRecordToFile = false;
                bPaused = true; //Pause So user Knows Saving Is disabled
            }
            if (!saveFoodTracks(vfishmodels,vfoodmodels,outfooddatafile,frameNumberString))
                bRecordToFile = false;
        }

        checkPauseRun(&window_main,keyboard,nFrame);


    } //main While loop
    //delete capture object
    capture.release();



    std::clog << gTimer.elapsed()/60000.0 << "[Progress] Exiting video processing loop <<<" <<std::endl;
    //Dont Forget to Reset startFrameCount = 1 So Next Video Starts from Beginning
    //stopFrame       = 0;//No Stopping on NExt Video
    //Close File
    closeDataFile(outdatafile);
    closeDataFile(outfooddatafile);

    return nFrame; //Return Number of Last Frame Processed
}




//Operator for Priority Ordering
bool operator<(const fishModel& a, const fishModel& b)
{
  return a.templateScore > b.templateScore; //Max Heap
}


///
/// \brief UpdateFishModels starting from Blob Info do the processing steps to update FishModels for this frame,
/// \param maskedImg_gray
/// \param vfishmodels
/// \param fishblobs
/// \param nFrame
/// \param frameOut
///\todo
///
void UpdateFishModels(const cv::Mat& maskedImg_gray,fishModels& vfishmodels,zftblobs& fishblobs,unsigned int nFrame,cv::Mat& frameOut){

    qfishModels qfishrank;

    fishModel* pfish = NULL;

    fishModels::iterator ft;
    bool bModelFound;

     // Look through Blobs find Respective fish model attached or Create New Fish Model if missing
    for (zftblobs::iterator it = fishblobs.begin(); it!=fishblobs.end(); ++it)
    {
        zftblob* fishblob = &(*it);

        cv::Point ptbcentre = fishblob->pt; //Start As First Guess - This is updated When TemplMatching
        cv::Point ptSearch; //Where To Centre The Template Matching Searcrh
        int bestAngle;
        double  maxMatchScore = 0.0;
        bModelFound = false;
        int iTemplRow = iLastKnownGoodTemplateRow; //Starting Search Point For Template
        int iTemplCol = 0;

        //Check Through Models And Find The Closest Fish To This FishBlob
        /// Note We do Template Matching On Previous Fish Coordinates First , (Not On Wobbly Blobs Coordinates)
        /// If No FishModel Is Matched with this Blob, then We Follow Up to Check the template score Of the Blob, before Creating A new Fish Model
        for ( ft  = vfishmodels.begin(); ft!=vfishmodels.end(); ++ft)
        {
            if (bModelFound) //Speed Up - Single Fish Mode - So Skip Others If Model Found
                continue;

             pfish = ft->second;
             assert(pfish);
             ///Does this Blob Belong To A Known Fish Model?
             //Check Overlap Of This Model With The Blob - And Whether The Image of this Blob contains something That looks like a fish
             if (pfish->zfishBlob.overlap(pfish->zfishBlob,*fishblob) > 0 )
             {
                 //If Yes then assign the fish with the overlapping blob the template Match Score
                bModelFound = true;
                ptSearch = pfish->ptRotCentre; //gptHead//((cv::Point)fishblob->pt-gptHead)/3+gptHead;
                iTemplRow = pfish->idxTemplateRow;
                iTemplCol = pfish->idxTemplateCol;
                maxMatchScore = doTemplateMatchAroundPoint(maskedImg_gray,ptSearch,iTemplRow,iTemplCol,bestAngle,ptbcentre,frameOut);

                //Failed? Try the blob Head (From Enhance Mask) Detected position
                if ( maxMatchScore < gTemplateMatchThreshold)
                {
                  ptSearch = ((cv::Point)fishblob->pt-gptHead)/3+gptHead;
                  maxMatchScore = doTemplateMatchAroundPoint(maskedImg_gray,ptSearch,iTemplRow,iTemplCol,bestAngle,ptbcentre,frameOut);
                }
                pfish->templateScore = maxMatchScore;
                //pfish->idxTemplateRow = iTemplRow; pfish->idxTemplateCol = iTemplCol;

                 if ( maxMatchScore >= gTemplateMatchThreshold)
                 {
                     //Some existing Fish Can be associated with this Blob - As it Overlaps from previous frame
                    ///Update Model State
                    // But not While it Is manually updating/ Modifying Bounding Box (Flags Are set in Mainwindow)
                    if (!bStoreThisTemplate && !bDraggingTemplateCentre) //Skip Updating Bound If this round we are saving The Updated Boundary
                    {
                        pfish->updateState(fishblob,maxMatchScore,bestAngle+iFishAngleOffset,ptbcentre,nFrame,gFishTailSpineSegmentLength,iTemplRow,iTemplCol);
                    }
                    else
                    { //Rotate Template Box - Since this cannot be done Manually
                        pfish->bearingAngle   = (bestAngle+iFishAngleOffset);
                        pfish->bearingRads   =  (bestAngle+iFishAngleOffset)*CV_PI/180.0;
                    }

                 }
                 else //Below Thres Match Score
                 {
                       //Overide If We cant find that fish anymore/ Search from the start of the row across all angles
                       if (pfish->inactiveFrames > gcMaxFishModelInactiveFrames)
                           iFishAngleOffset = 0;
                         qDebug() << nFrame << " Guessing next TemplCol:" << iFishAngleOffset;
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
            //Check Template Match Score
            ptSearch = fishblob->pt;
            // Suitable template has not been found yet for this model, starting search Point For Template can be were we left off, as a suitable
             iTemplRow = iLastKnownGoodTemplateRow ;
             iTemplCol = 0;
            pwindow_main->LogEvent("No Fish model found for blob");
            cv::circle(frameOut,ptSearch,3,CV_RGB(15,15,250),1); //Mark Where Search Is Done
            maxMatchScore = doTemplateMatchAroundPoint(maskedImg_gray,ptSearch,iTemplRow,iTemplCol,bestAngle,ptbcentre,frameOut);
            //If New Blob Looks Like A Fish, Then Make  A New Model
            if (maxMatchScore > gTemplateMatchThreshold*0.90)
            {
                //Make new fish Model
               fishModel* fish= new fishModel(*fishblob,bestAngle,ptbcentre);
               fish->ID = ++gi_MaxFishID;
               fish->idxTemplateRow = iTemplRow;
               fish->idxTemplateCol = iTemplCol;
               fish->updateState(fishblob,maxMatchScore,bestAngle,ptbcentre,nFrame,gFishTailSpineSegmentLength,iTemplRow,iTemplCol);

               vfishmodels.insert(IDFishModel(fish->ID,fish));
               qfishrank.push(fish); //Add To Priority Queue
               std::stringstream strmsg;
               strmsg << " New fishmodel: " << fish->ID << " with Template Score :" << fish->templateScore;
               //std::clog << nFrame << strmsg.str() << std::endl;
               pwindow_main->LogEvent(QString::fromStdString(strmsg.str()));
            }

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
                if (!iroi.contains(pfishBest->ptRotCentre,pfishBest->bodyRotBound.size.width+pfishBest->bodyRotBound.size.height ))
                {
                    qfishrank.pop();
                    pfishBest = 0;
                }
             }
   }//Search For Best Model

   if (pfishBest)
   {
        //qfishrank.pop();//Remove From Priority Queue Rank
        maxTemplateScore = pfishBest->templateScore;
        pfishBest->inactiveFrames   = 0; //Reset Counter
    }

   /// Delete All FishModels EXCEPT the best Match - Assume 1 Fish In scene / Always Retain 1 Model
    ft = vfishmodels.begin();
    while(ft != vfishmodels.end() ) //&& vfishmodels.size() > 1
    {
        pfish = ft->second;
        // If this is not the Best One
        if (pfishBest != pfish ) //&& pfishBest != 0
        {
            //Check Ranking Is OK, as long off course that a fishTemplate Has Been Found On This Round -
            //OtherWise Delete The model?
            //Assertion Fails When Old Model Goes Out Of scene and video Is retracked
            //assert(pfish->templateScore < maxTemplateScore || maxTemplateScore == 0);
            //If We found one then Delete the other instances waiting for a match - 1 Fish Tracker
            if (bModelFound || pfish->inactiveFrames > gcMaxFishModelInactiveFrames) //Check If it Timed Out / Then Delete
            {
                std::clog << gTimer.elapsed()/60000 << " " << nFrame << "# Deleted fishmodel: " << pfish->ID << "Inactive:"<< pfish->inactiveFrames << " Low Template Score :" << pfish->templateScore << " when Best is :"<< maxTemplateScore << std::endl;
                ft = vfishmodels.erase(ft);
                delete(pfish);
                continue;
            }else
            // Fish Model Has not Been Matched / Increment Time This Model Has Not Been Active
                 pfish->inactiveFrames ++; //Increment Time This Model Has Not Been Active

        }
        ++ft; //Increment Iterator
    } //Loop To Delete Other FishModels



} //UpdateFishModels //


//        cv::Point pBound2 = cv::Point(max(0,min(maskedImg_gray.cols,centroid.x+gFishBoundBoxSize+2)), max(0,min(maskedImg_gray.rows,centroid.y+gFishBoundBoxSize+2)));



/// Process Optic Flow of defined food model positions
/// Uses Lukas Kanard Method to get the estimated new position of Prey Particles
///
int processFoodOpticFlow(const cv::Mat frame_grey,const cv::Mat frame_grey_prev,foodModels& vfoodmodels,unsigned int nFrame,zftblobs& vPreyKeypoints_next )
{
    int retCount = 0;
   std::vector<cv::Point2f> vptPrey_current;
   std::vector<cv::Point2f> vptPrey_next;


    zftblobs vPreyKeypoints_current;
    zftblobs vPreyKeypoints_ret;
   std::vector<uchar> voutStatus;
   // L1 distance between patches around the original and a moved point, divided by number of pixels in a window, is used as a error measure.
   std::vector<float>    voutError;

   foodModel* pfood = NULL;
   foodModels::iterator ft;

    //Fill POint Vector From foodmodel vector
   for ( ft  = vfoodmodels.begin(); ft!=vfoodmodels.end(); ++ft)
   {
       pfood = ft->second;
       cv::KeyPoint kptFood(pfood->zTrack.centroid,pfood->zfoodblob.size);
       vPreyKeypoints_current.push_back(kptFood  );
   }

    cv::KeyPoint::convert(vPreyKeypoints_current,vptPrey_current);

    //Calc Optic Flow for each food item
    if (vptPrey_current.size() > 0 && !frame_grey_prev.empty())
        cv::calcOpticalFlowPyrLK(frame_grey_prev,frame_grey,vptPrey_current,vptPrey_next,voutStatus,voutError,cv::Size(31,31),2);

    cv::KeyPoint::convert(vptPrey_next,vPreyKeypoints_next);

    //update food item Location
        //Loop through points
    for (int i=0;i<(int)vPreyKeypoints_ret.size();i++)
    {
        if (!voutStatus.at(i))
            continue; //ignore bad point
        vPreyKeypoints_next.push_back(vPreyKeypoints_ret.at(i)); //fwd the good ones
//        // find respective food model, update state
//        vfoodmodels[i]->zTrack.centroid = vPreyKeypoints_next.at(i);
//        vfoodmodels[i]->zfoodblob.pt = vPreyKeypoints_next.at(i);
//        vfoodmodels[i]->updateState(&vfoodmodels[i]->zfoodblob,0,
//                                    vPreyKeypoints_next.at(i),
//                                    nFrame,vfoodmodels[i]->blobMatchScore,
//                                    vfoodmodels[i]->blobRadius);
//        retCount++;
    } //Check if Error
//
return retCount;
}


///
/// \brief The foodBlobMatch struct used to hold the score when matching blob to food item.
/// It allows score vector to be sorted
///
struct foodBlobMatch
{
    int score;
    foodModel* pFoodObject;
    zfdblob*   pFoodBlob;

    foodBlobMatch(foodModel* pFood,  zfdblob* pBlob,int inscore) : pFoodObject(pFood), pFoodBlob(pBlob),score(inscore) {}

    bool operator < (const foodBlobMatch& pair) const
    {
        return (score > pair.score);
    }
};


///
/// \brief UpdateFoodModels A score based assignment of blob to model that uses overlap, distance and blob size- comparison
/// between foodobject laststate and new frame blobs to update food state on the most likely match to blob
/// Can be converted to statistical model of assignment
/// \param maskedImg_gray
/// \param vfoodmodels
/// \param foodblobs
/// \param nFrame
/// param frameOut (removed) - no drawing should happen here
/// \todo Add calcOpticalFlowPyrLK Lucas-Kanard Optic Flow Measurment to estimate food displacement
void UpdateFoodModels(const cv::Mat& maskedImg_gray,foodModels& vfoodmodels,zfdblobs& foodblobs,unsigned int nFrame,bool bAllowNew=true)
{
    qfoodModels qfoodrank;
    foodModel* pfood = NULL;

    foodModels::iterator ft;
    //Make Triplet of Score, and foodModel, zfdBlob
    std::vector<foodBlobMatch> vPairScores;
    zfdblobs vfoodblobs_spare; //A copy of foodblobs from which we delete the blobs that we match to food objects
    //static const int nfoodBlobs  = G_MAX_FOOD_BLOB_LIMIT;
    //static const int nfoodObjects = G_MAX_FOOD_BLOB_LIMIT;
    //Make score Array blobsXFoodModels - Fixed Size
    //static double dFoodScoreToBlob[nfoodBlobs][nfoodObjects];//Vector of \nabla d for error functions
    //memset(dFoodScoreToBlob,0.0,nfoodBlobs*nfoodObjects*sizeof(double));

    int bIdx,fIdx;

        /// Assign Blobs To Food Models //
     // Look through KeyPoints/Blobs find Respective food model attached or Create New Food Model if missing
    bIdx = 0;
    for (zfdblobs::iterator it = foodblobs.begin(); it!=foodblobs.end(); ++it)
    {
        zfdblob* foodblob  = &(*it);
        cv::Point2f ptblobCentroid = foodblob->pt;

        fIdx = 0;

        //Score Distance To Each Food
        for ( ft  = vfoodmodels.begin(); ft!=vfoodmodels.end(); ++ft)
        {
             pfood = ft->second;
            pfood->blobMatchScore = 0 ;
            cv::Point2f ptfoodCentroid = pfood->zTrack.centroid;

            //Initial score Is high and drops with penalties
            int iMatchScore = gMaxClusterRadiusFoodToBlob*2; //Reset to 5So We Can Rank this Match

            float overlap = pfood->zfoodblob.overlap(pfood->zfoodblob,*foodblob);
            float fbdist = norm(ptfoodCentroid-ptblobCentroid);

            ///TODO Can add Directional comparison - motion accross frames should be similar
            /// Also a gaussian model of speed, size, direction - tuned to each foodobject could be nice here

            //Penalize Distance in px moved since last seen

            iMatchScore -= exp(2* fbdist/ (nFrame - pfood->nLastUpdateFrame + 1) ) ;

            //Penalize Size Mismatch
            iMatchScore -= 2.0*abs(pfood->zfoodblob.size - foodblob->size);

            //Bonus For Overlap
            if (overlap > 0.0)
                iMatchScore +=(int)(10.0*overlap);

            //Only add Pairs with possible Scores - (Set By the initial iMatchScore Value)
            if (iMatchScore > 0)
            {
                foodBlobMatch pBlobScore(pfood,foodblob,iMatchScore);

                // qDebug() << "Invalid blob-prey pair detected";
                //float d = cv::norm(pfood->zTrack.centroid - foodblob->pt);
                //assert(d > gMaxClusterRadiusFoodToBlob)
                vPairScores.push_back(pBlobScore); //Append Pair Score to vector
            }

            fIdx++;
         } // Loop Through Food Models

        bIdx++;


        //Safety Net
       // if (vPairScores.size() > 1000)
        //    break;
    } // Loop Through BLOBS


   ///
   /// Run Through Possible Food-KeyPoint Pairs and update models of best matches  //
   // Sort Score Vector so we hit best food-blob match first //
   std::sort(vPairScores.begin(),vPairScores.end() );
   int iMatchCount =0;
   zfdblob foodblobMatched;

   //Make A copy of the Blob vector So we do not Mess up the Pair.blob pointers
   vfoodblobs_spare = foodblobs;
   // Start from Top Score And Assign pairs by rank //
   std::vector<foodBlobMatch>::iterator it = vPairScores.begin();
   while( it!=vPairScores.end())
   {
        foodBlobMatch pair = (*it);
        /// Check if Blob Has Already been taken out And Matched //
        bool blobAvailable = false; //Flag If this Blob Has been Preivously Matched to A food (Ie Used)
        //bool blobConsumed = false; //Flag If this Blob Has Now been and should be deleted

        foodblobMatched = *pair.pFoodBlob;
        //Find matching keypoint/ and remove from list of available blobs for matching
        blobAvailable = checkKeypointExistsAndRemove(vfoodblobs_spare,foodblobMatched);

        //Candidate Pair, directs (Points) to Blob Which is not Available anymore /already matched blob
        if (!blobAvailable)
        {   ++it;
            continue; //Check next Pair
        }

        //Skip Paired Up FoodObject Items
        if (blobAvailable && (pair.pFoodObject->nLastUpdateFrame - nFrame) > 0 )
        {
            // Update Food Item Based oN best Blob Match
            pair.pFoodObject->activeFrames ++; //Increase Count Of Consecutive Active Frames
            pair.pFoodObject->updateState(foodblobMatched,0, foodblobMatched.pt,nFrame, pair.score,foodblobMatched.size);
            iMatchCount++;
           // blobConsumed = true;
        }
       ++it;
   }//Loop Through Pairs And MAtch Food To Blob

//qDebug() << "Matched " << iMatchCount;
    /// Unmatched Keypoints/blobs handling -create new Food Models ///
    // Check For Unmatched Blobs And Make (A) New Food Item Each Time this function is called //
    //If There Are more Blobs Remaining / New Objects Allowed (Currently OFF for OpticFlow, and we are below the Count Limit
    if (bAllowNew && (vfoodblobs_spare.size() > 1) && (vfoodmodels.size() < gi_FoodModelNumberLimit) )
    {
        zfdblob* foodblob= &vfoodblobs_spare[0]; //Take the 1st one Available And Make A food Model For it
        pfood = pwindow_main->getFoodItemAtLocation(foodblob->pt); //Check for dublicated
        if (!pfood) //IF no dublicate Item there, then make new
        {
            pfood = new foodModel(*foodblob ,++gi_MaxFoodID);
            pfood->blobMatchScore = 0;
            vfoodmodels.insert(IDFoodModel(pfood->ID,pfood));
            //_DEBUG
            //std::stringstream strmsg;
            //strmsg << "# New foodmodel: " << pfood->ID << " N:" << vfoodmodels.size();
            //std::clog << nFrame << strmsg.str() << std::endl;
            pfood->updateState(*foodblob,0,foodblob->pt,nFrame,500,foodblob->size);
        }else{
            //qDebug() << "Dublicate prey location for ID: " << pfood->ID ;
        }
    } //Make New Food Model If Allowed


    /// Delete All Inactive Food Models ///
    ft = vfoodmodels.begin();
    while(ft != vfoodmodels.end() ) //&& vfishmodels.size() > 1
    {
        pfood = ft->second;
        //Delete If Not Active for Long Enough between inactive periods / Track Unstable
         if (pfood->isUnused())
        {
            ft = vfoodmodels.erase(ft);
            //_DEBUG
            //std::clog << nFrame << "# Delete foodmodel: " << pfood->ID << " Inactive:" << pfood->inactiveFrames <<  " N:" << vfoodmodels.size() << std::endl;
            delete(pfood);

            continue;
        }
        else //INcrease Inactive Frame Count For this Food Model
        {//If this Model Has not Been Used Here
            if (pfood->nLastUpdateFrame-nFrame > 1)
            {
                //pfood->activeFrames = 0; //Reset Count Of Consecutive Active Frames
                pfood->inactiveFrames ++; //Increment Time This Model Has Not Been Active
            }
        }

        ++ft; //Increment Iterator
    } //Loop To Delete Inactive FoodModels

} //UpdateFoodModels //


bool checkKeypointExistsAndRemove(zfdblobs& vfoodblobs_spare,zfdblob& FoodBlob)
{
    bool blobAvailable = false;
    //zfdblob foodblobMatched;
    ///Find Paired Blob From Available Blobs - If There then Still Available
    /// this is done so we Can detect the unmatched Ones
    zfdblobs::iterator ft = vfoodblobs_spare.begin();
    while(ft != vfoodblobs_spare.end() )
    {
        zfdblob* foodblob = &(*ft);
        //is this blob/keypoint the same / ie still available in spare list?
        if ( cv::norm(foodblob->pt - FoodBlob.pt ) < 0.3 )
        {
            //Keep a copy before deleting
            //foodblobMatched = *foodblob;
            blobAvailable = true;
            //Deleting Blobs Invalidates the next Pointer Pairs! So do not alter the vector Just yet!
            //if (blobConsumed) //Erase Matched Blob From Available List - This can happen when Paused
            ft = vfoodblobs_spare.erase(ft);

            break; //Found Blob -> Exit Loop
        }
        else
            ++ft;


    }//Loop And Find Used Blob

    return blobAvailable;
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
        pwindow_main->LogEvent("[info] User Terminated Tracker- Bye!");
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
            iFishAngleOffset = 0;
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
        iFishAngleOffset--;
        iFishAngleOffset = max(-180,iFishAngleOffset);
        pwindow_main->LogEvent(QString("User Rotated Template:")+QString::number(iFishAngleOffset)  );
    }

    if ((char)keyboard == ']') //Rotate Template ClockWise
    {
        iFishAngleOffset++;
        iFishAngleOffset = min(180, iFishAngleOffset);
        pwindow_main->LogEvent(QString("User Rotated Template:")+QString::number(iFishAngleOffset)  );
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
                 fish->c_spineSegL = gFishTailSpineSegmentLength;
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


    if ((char)keyboard == 'w')
    {
      bRecordToFile = !bRecordToFile; //Main Loop Will handle this
      if (bRecordToFile)
      {
        pwindow_main->LogEvent(QString(">> [DISABLED] Recording Tracks ON - New File <<"));
        ///Code Moved TO MainWindow GUI
      }
      else
        pwindow_main->LogEvent(QString("<< [DISABLED] Recording Tracks OFF >>"));
    }

    /// Manual Prey addition mode - a left click adds new prey
    if ((char)keyboard == 'P')
    {

        bAddPreyManually = !bAddPreyManually;

        if (bAddPreyManually)
            pwindow_main->LogEvent(QString(">> Manual Prey Adding ON <<"));
        else
            pwindow_main->LogEvent(QString("<< Manual Prey Adding OFF >>"));
    }

    /// Measure Distance In straight Line Mode
    if ((char)keyboard == 'M')
    {
       bMeasure2pDistance= !bMeasure2pDistance;
       if (bMeasure2pDistance)
       {
           pwindow_main->LogEvent(QString(">> Manual Distance Measurement ON <<"));
           pwindow_main->SetTrackerState(7);
           bPaused = true;
       }
       else
           pwindow_main->LogEvent(QString("<< Manual Distance Measurement  OFF >>"));
           //pwindow_main->ui->statusBar->showMessage(("Adjust fish detection template position"));
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
        iFishAngleOffset = 0;
    }

    if ((char)keyboard == 'x')
    {
        gUserReward = -1000;
        pwindow_main->LogEvent("User -ve Reward ");
    }
    if ((char)keyboard == 'a')
    {
        gUserReward = +1000;
        pwindow_main->LogEvent("User +ve Reward ");
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

    dirToSave.append("/pics/" + fileVidCoreName + "/");
    //QString imageToSave =  fileVidCoreName + "_" + frameNumberString + ".png";
    QString imageToSave = frameNumberString + ".png";
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

            if (iroi.contains(kp.pt,gszTemplateImg.width ))
                     ptFishblobs.push_back(kp);

            //int maskVal=(int)gframeMask.at<uchar>(kp.pt);
            //if (maskVal > 0)
             //keypoints_in_mask.push_back(kp);
        }
    }


    // Draw detected blobs as red circles.
    // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures the size of the circle corresponds to the size of blob
    //frame.copyTo(frameOut,maskimg); //mask Source Image
    //cv::drawKeypoints( frameOut, ptFishblobs, frameOut, cv::Scalar(250,20,20), cv::DrawMatchesFlags::DEFAULT ); //cv::DrawMatchesFlags::DRAW_RICH_KEYPOINTS


    detector->clear();

}


/// Updated Blob Processing
/// \brief processFoodBlobs Finds blobs that belong to rotifers
/// \param frame
/// \param maskimg
/// \param frameOut //Output Image With FishBlob Rendered
/// \param ptFoodblobs opencv keypoints vector of the Fish
/// \return
/// \note Draws Blue circle around food blob, with relative size
///
int processFoodBlobs(const cv::Mat& frame_grey,const cv::Mat& maskimg,cv::Mat& frameOut,std::vector<cv::KeyPoint>& ptFoodblobs)
{

    cv::Mat frameMasked;
    //Thresholds Set By Gui
    //g_SegFoodThesMax = g_Segthresh*1.25; //Set to current Seg Thresh
    //g_SegFoodThesMin = g_Segthresh*0.80;

    if (!maskimg.empty())
        frame_grey.copyTo(frameMasked,maskimg); // Do not Apply Mask
    else
        frame_grey.copyTo(frameMasked);

    std::vector<cv::KeyPoint> keypoints;


    //std::vector<cv::KeyPoint> keypoints_in_ROI;
    cv::SimpleBlobDetector::Params params;

    //a circle has a circularity of 1,
    //circularity of a square is 0.785, and so on.

    params.filterByCircularity  = true;
    params.minCircularity       = 0.60;
    params.maxCircularity       = 1.0;

    params.filterByColor        = false;
    params.filterByConvexity    = false;


    params.maxThreshold = g_SegFoodThesMax; //Use this Scanning to detect smaller Food Items
    params.minThreshold = g_SegFoodThesMin;
    params.thresholdStep = 4;

    // Filter by Area.
    params.filterByArea = true;
    params.minArea = 0;
    params.maxArea = gthres_maxfoodblobarea;

    /////An inertia ratio of 0 will yield elongated blobs (closer to lines)
    ///  and an inertia ratio of 1 will yield blobs where the area is more concentrated toward the center (closer to circles).
    params.filterByInertia      = false;
    params.maxInertiaRatio      = 1.0;
    params.minInertiaRatio      = 0.1;

    params.minDistBetweenBlobs = gMaxClusterRadiusFoodToBlob/2;

    //params.filterByInertia = true;

    // Set up the detector with default parameters.
    cv::Ptr<cv::SimpleBlobDetector> detector = cv::SimpleBlobDetector::create(params);

    assert(frameMasked.depth() == CV_8U);
    detector->detect( frameMasked, keypoints,maskimg); //frameMask


    //Mask Is Ignored so Custom Solution Required
    //for (cv::KeyPoint &kp : keypoints)

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
    if (bDrawFoodBlob)
        cv::drawKeypoints( frameOut, ptFoodblobs, frameOut, cv::Scalar(0,120,200), cv::DrawMatchesFlags::DEFAULT );


    detector->clear();

    return (int)ptFoodblobs.size();

}

//int saveTrackedBlobs(cvb::CvBlobs& blobs,QString filename,std::string frameNumber,ltROI& roi)
//{
//    int cnt = 0;
//    int Vcnt = 1;
//    bool bNewFileFlag = true;

//    //Loop Over ROI
//    Vcnt++; //Vial Count
//    cnt = 0;

//    QFile data(filename);
//    if (data.exists())
//        bNewFileFlag = false;

//    if(data.open(QFile::WriteOnly |QFile::Append))
//    {
//        QTextStream output(&data);
//        if (bNewFileFlag)
//             output << "frameN,SerialN,BlobLabel,Centroid_X,Centroid_Y,Area\n" ;

//        //Loop Over Blobs
//        for (cvb::CvBlobs::const_iterator it=blobs.begin(); it!=blobs.end(); ++it)
//        {

//            cvb::CvBlob* cvB = it->second;
//            cv::Point pnt;
//            pnt.x = cvB->centroid.x;
//            pnt.y = cvB->centroid.y;

//            cnt++;

//            if (roi.contains(pnt))
//                //Printing the position information
//                output << frameNumber.c_str() << "," << cnt <<","<< cvB->label << "," << cvB->centroid.x <<","<< cvB->centroid.y  <<","<< cvB->area  <<"\n";
//          }


//       data.close();

//      }


//    return cnt;
//}

////Saves the total Number of Counted Blobs and Tracks only
//int saveTrackedBlobsTotals(cvb::CvBlobs& blobs,cvb::CvTracks& tracks,QString filename,std::string frameNumber,ltROI& roi)
//{

//    bool bNewFileFlag = true;
//    int cnt = 0;
//    int Larvacnt = 0;
//    cnt++;
//    //cv::Rect iroi = (cv::Rect)(*it);

//    QFile data(filename);
//    if (data.exists())
//        bNewFileFlag = false;

//    if(data.open(QFile::WriteOnly |QFile::Append))
//    {

//        int blobCount = 0;
//        int trackCount = 0;

//        //Count Blobs in ROI
//        for (cvb::CvBlobs::const_iterator it = blobs.begin(); it!=blobs.end(); ++it)
//        {
//            cvb::CvBlob* blob = it->second;
//            cv::Point pnt;
//            pnt.x = blob->centroid.x;
//            pnt.y = blob->centroid.y;

//            if (roi.contains(pnt))
//                blobCount++;
//        }

//        //Count Tracks in ROI
//        for (cvb::CvTracks::const_iterator it = tracks.begin(); it!=tracks.end(); ++it)
//        {
//            cvb::CvTrack* track = it->second;
//            cv::Point pnt;
//            pnt.x = track->centroid.x;
//            pnt.y = track->centroid.y;

//            if (roi.contains(pnt))
//                trackCount++;
//        }


//        QTextStream output(&data);
//        if (bNewFileFlag)
//             output << "frameN,blobN,TracksN \n";

//        output << frameNumber.c_str() << "," << blobCount << "," << trackCount <<"\n";
//        Larvacnt +=blobCount;
//        data.close();
//    }


//    return Larvacnt;
//}


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

///
/// \brief resetDataRecording Clean the Output File, And Starts Over -
/// Triggered when Recording Is toggled on - such that a fresh file is created Each Time
/// \param strpostfix / Either food or tracks, added to the file name with a sequential Number
/// \return  True if file opened Succesfully
///
bool resetDataRecording(QFile& outdatafile,QString strpostfix)
{
    closeDataFile(outdatafile); //
    //removeDataFile(outdatafile);
    //extract Post Fix
    QFileInfo fileInfFish(outdatafile);

    if ( !openDataFile(fileInfFish.absoluteDir().absolutePath(),fileInfFish.completeBaseName(),outdatafile,strpostfix) )
    {
        pwindow_main->LogEvent(QString("[Error] Opening Data " + strpostfix +" Tracks File"));
        return false;
    }

    return true;

}

void writeFishDataCSVHeader(QFile& data)
{

    /// Write Header //
    QTextStream output(&data);
    output << "frameN \t ROI \t fishID \t AngleDeg \t Centroid_X \t Centroid_Y \t EyeLDeg \t EyeRDeg \t ThetaSpine_0 \t ";
    for (int i=1;i<gFishTailSpineSegmentCount;i++)
        output <<  "DThetaSpine_" << i << "\t";

    output << " tailSegmentLength";
    output << "\t templateScore";
    output << "\t lastTailFitError";
    output << "\t lEyeFitScore";
    output << "\t rEyeFitScore";
    output << "\t nFailedEyeDetectionCount";
    output << "\t RotiferCount \n";

}


void writeFoodDataCSVHeader(QFile& data)
{
    /// Write Header //
    QTextStream output(&data);
    output << "FrameN \t ROI \t FoodID \t Centroid_X \t Centroid_Y \t Radius \t InactiveFrames \n";

}


bool openDataFile(QString filepathCSV,QString filenameVid,QFile& data,QString strpostfix)
{
    int Vcnt = 1;
    bool newFile = false;
    //Make ROI dependent File Name
    QFileInfo fiVid(filenameVid);
    QFileInfo fiOut(filepathCSV+"/") ;
    QString fileVidCoreName = fiVid.completeBaseName();
    QString dirOutPath = fiOut.absolutePath() + "/"; //filenameCSV.left(filenameCSV.lastIndexOf("/")); //Get Output Directory

    //strpostfix = strpostfix + "_%d.csv";


    //char buff[50];
    //sprintf(buff,strpostfix.toStdString(),Vcnt);
    //dirOutPath.append(fileVidCoreName); //Append Vid Filename To Directory
    //dirOutPath.append(buff); //Append extension track and ROI number
    if (fileVidCoreName.contains(strpostfix,Qt::CaseSensitive))
    {
        fileVidCoreName = fileVidCoreName.left(fileVidCoreName.lastIndexOf("_"));
        dirOutPath = dirOutPath + fileVidCoreName+ "_" + QString::number(Vcnt) +  ".csv";
    }
    else
        dirOutPath = dirOutPath + fileVidCoreName + strpostfix + "_" + QString::number(Vcnt) + ".csv";

    data.setFileName(dirOutPath);
    //Make Sure We do not Overwrite existing Data Files
    while (!newFile)
    {
        if (!data.exists() || data.isOpen()) //Write HEader
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
            //Filename Is Like AutoSet_12-10-17_WTNotFedRoti_154_002_tracks_1.csv
                //Increase Seq Number And Reconstruct Name
                Vcnt++;
                // If postfix (track / food) already there, then just add new number
                if (fileVidCoreName.contains(strpostfix,Qt::CaseSensitive))
                    dirOutPath = fiOut.absolutePath() + "/" + fileVidCoreName + "_" + QString::number(Vcnt) + ".csv";
                else
                    dirOutPath = fiOut.absolutePath() + "/" + fileVidCoreName + strpostfix + "_" + QString::number(Vcnt) + ".csv";



                data.setFileName(dirOutPath);
                //data.open(QFile::WriteOnly)

            }
         }
    }
    if (!data.open(QFile::WriteOnly |QFile::Append))
    {
        std::cerr << "Could not open output file : " << data.fileName().toStdString() << std::endl;
        return false;
    }else {
        //New File
        if (!bBlindSourceTracking)
        std::clog << "Opened file " << dirOutPath.toStdString() << " for data logging." << std::endl;

        //output.flush();

    }

    return true;
}




void closeDataFile(QFile& data)
{
    data.close();
    if (bBlindSourceTracking)
        std::clog << gTimer.elapsed()/60000 << " Closed Output File " << std::endl;
    else
        std::clog << gTimer.elapsed()/60000 << " Closed Output File " << data.fileName().toStdString() << std::endl;
}

void removeDataFile(QFile& data)
{
    if (bBlindSourceTracking)
        std::clog << gTimer.elapsed()/60000 << "[Warning] Deleting Output File " << std::endl;
    else
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
int saveTracks(fishModels& vfish,foodModels& vfood,QFile& fishdata,QString frameNumber)
{
    int cnt;
    int Vcnt = 0;


    //Make ROI dependent File Name
    if (!fishdata.exists())
    {
        std::cerr << "Fish Track File Is Missing" << std::endl;
        pwindow_main->LogEvent("[Error] Fish Track File Is Missing");
        return 0;
    }

    //Loop Over ROI
    for (ltROIlist::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
    {
        cnt = 1;
        Vcnt++;
        ltROI iroi = (ltROI)(*it);
        //Make ROI dependent File Name

        QTextStream output(&fishdata);

        //Save Tracks In ROI
        for (fishModels::iterator it=vfish.begin(); it!=vfish.end(); ++it)
        {
            cnt++;
            fishModel* pfish = it->second;
            //cvb::CvLabel cvL = it->first;

            if (iroi.contains(pfish->ptRotCentre))
            {
                //Printing the position information +
                //+ lifetime; ///< Indicates how much frames the object has been in scene.
                //+active; ///< Indicates number of frames that has been active from last inactive period.
                //+ inactive; ///< Indicates number of frames that has been missing.
                output << frameNumber << "\t" << Vcnt  << "\t" << (*pfish);
                output << "\t" << foodModel::getActiveFoodCount(vfood) << "\n";
            }
            //Empty Memory Once Logged
            pfish->zTrack.pointStack.clear();
            pfish->zTrack.pointStack.shrink_to_fit(); //Requires this Call of C++ otherwise It Doesnt clear
         }//For eAch Fish

        //Regular Timed entry - in the absence of fish
         //If there is are no fish Then Add a regular Entry denoting the number of prey
        if (bTrackFood && vfish.size() == 0 && (frameNumber.toUInt()%gFoodReportInterval == 0 || frameNumber.toUInt()==1))
        {
            //make Null Fish
            fishModel* pNullfish = new fishModel();
            pNullfish->ID          = 0;
            pNullfish->resetSpine();

            output << frameNumber << "\t" << Vcnt  << "\t" << (*pNullfish);
            output << "\t" << foodModel::getActiveFoodCount(vfood) << "\n";
            delete pNullfish;
        }

   } //Loop ROI

 return cnt;
} //saveTracks

int saveFoodTracks(fishModels& vfish,foodModels& vfood,QFile& fooddata,QString frameNumber)
{

    //Make ROI dependent File Name
    if (!fooddata.exists())
    {
        std::cerr << "Prey Model File Is Missing" << std::endl;
        pwindow_main->LogEvent("[Error] Prey File Is Missing");
        return 0;
    }

    QTextStream output(&fooddata);

    foodModels::iterator ft = vfoodmodels.begin();
    while (ft != vfoodmodels.end())
    //for (int i =0;i<v.size();i++)
    {
        foodModel* pfood = ft->second;

         if (pfood->isTargeted) //Only Log The Marked Food
         {
            output << frameNumber << "\t" << pfood->ROIID << "\t" << pfood->ID << "\t" << pfood->zTrack << "\t"  << pfood->blobRadius << "\t"  << pfood->inactiveFrames << "\n";

            pfood->zTrack.pointStack.clear();
            pfood->zTrack.pointStack.shrink_to_fit();
         }
    ++ft;
    }





    return 1;
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

/// Draws all paired points found in the vMeasureLines vector - where user defined points are stored by the user when
/// in measuremode .
void drawUserDefinedPoints(cv::Mat& frame)
{
    for (std::vector<pointPair>::iterator it = vMeasureLines.begin(); it != vMeasureLines.end(); ++it)
    {
        pointPair pPoints = (pointPair)(*it);
        cv::circle(frame,pPoints.first,1,cv::Scalar(0,200,60),1);
        cv::circle(frame,pPoints.second,1,cv::Scalar(20,60,220),1);
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


/////
///// \brief findMatchingContour Looks for the inner contour in a 2 level hierarchy that matches the point coords
///// \param contours source array in which to search
///// \param hierarchy
///// \param pt - Position around which we are searching
///// \param level - The required hierarchy level description of the contour being searched for
///// \param matchhull approx shape we are looking for
///// \param fittedEllipse Not Used - pointer to array of Rotated rect fitted ellipsoids
///// \return Index of *child*/Leaf contour closest to point
/////
//int findMatchingContour(std::vector<std::vector<cv::Point> >& contours,
//                              std::vector<cv::Vec4i>& hierarchy,
//                              cv::Point pt,
//                              int level,
//                              std::vector<cv::Point>* matchhull = nullptr,
//                              std::vector<cv::RotatedRect>* outfittedEllipse = nullptr)
//{
//    int idxContour           = -1;
//    bool bContourfound       = false;
//    int mindistToCentroid    = +10000; //Start Far
//    int distToCentroid       = +10000;
//    int matchContourDistance = 10000;

//    int tArea = 0;
//    int sArea = 0;

//    int tLength = 0;
//    int sLength = 0;

//    double dHudist = 0.0; //Shape Distance Hu moments distance measure from OpenCV

//    /// Render Only Countours that contain fish Blob centroid (Only Fish Countour)
//   ///Search Through Contours - Draw contours + hull results

//   ///Find Contour with Min Distance in shape and space -attach to closest contour
//   //In Not found Search Again By distance tp Full Contour
//       //Find Closest Contour
//       for( int i = 0; i< (int)contours.size(); i++ )
//       {

//          //Filter According to desired Level
//          if (level == 0) /////Only Process Parent Contours
//          {
//            if (hierarchy[i][3] != -1) // Need to have no parent
//               continue;
//            if (hierarchy[i][2] == -1)  // Need to have child
//                continue;
//            assert(hierarchy[hierarchy[i][2]][3] == i ); // check that the parent of the child is this contour i
//          }

//          if (level == 1) /////Only Process Child Contours
//          {
//              if (hierarchy[i][3] == -1) // Need to have a parent
//                  continue;
////                   //Parent should be root
////                   if (hierarchy[hierarchy[i][3]][3] != -1)
////                       continue;
//          }

//          if (level == 2) ////Needs to be top Level Contour
//          {
//              if (hierarchy[i][3] != -1) // No Parent Contour
//                  continue;
////                   //Parent should be root
////                   if (hierarchy[hierarchy[i][3]][3] != -1)
////                       continue;
//          }



//          //It returns positive (inside), negative (outside), or zero (on an edge)
//          //Invert Sign and then Rank From Smallest to largest distance
//          if (contours[i].size() > 0)
//            matchContourDistance = distToCentroid = -cv::pointPolygonTest(contours[i],pt,true);

//          //Measure Space Mod -Penalize Outside contour Hits - Convert Outside Contour Distances to X times further +ve (penalize)
//          //Make Distance alway positive
//          //matchContourDistance = (distToCentroid<0)?abs(distToCentroid)*20:distToCentroid;
//          // qDebug() << "-c" << i << " D:" <<  distToCentroid;


//          ///Match Shape -
//          /// \warning  If initial Shape Is not eye like this may be stuck into rejecting shapes
//          // If A shape is provided
//          //Find Contour Shape Similar to the one previously used for eye(ellipsoid)
//          if (matchhull != nullptr)
//          {
//              if (matchhull->size() > 5 && gOptimizeShapeMatching) //Only If Shape has been initialized/Given
//              {
//                   dHudist = cv::matchShapes(*matchhull,contours[i],CV_CONTOURS_MATCH_I2,0.0);
//                   matchContourDistance += dHudist*10.0; //Add Shape Distance onto / X Scale so it obtains relative importance
//                   // Now Check That distance is not too far otherwise reject shape
//                   //if (matchContourDistance > 1.0)
//                   //    continue; //Next Shape/Contour
//                   //qDebug() << "HuDist:" << dHudist*10.0;
//                 //Check Area

//                   tArea = cv::contourArea(*matchhull);
//                   sArea = cv::contourArea(contours[i]);

//                   tLength = cv::arcLength(*matchhull,true);
//                   sLength = cv::arcLength(contours[i],true);

//                   //Add Difference in Area to Score
//                   matchContourDistance += (abs(tArea - sArea));
//                  // qDebug() << "AreaDist:" << abs(tArea - sArea);

//                   matchContourDistance += abs(tLength - sLength);
//                   //qDebug() << "LengthDist:" << abs(tLength - sLength);

//              }
//          } // If Hull To Search For is provided


//          //Only Update if Spatial Distance is smaller but also moving from outside to inside of the shape
//          //-ve is outside - 0 on border -
//          //if(mindistToCentroid <= 0 && distToCentroid >= 0))
//          {
//               if (matchContourDistance < mindistToCentroid)
//               {
//                   //Otherwise Keep As blob Contour
//                   idxContour = i;
//                   mindistToCentroid = matchContourDistance;//New Min

//                   //qDebug() << "-----MinTD:"<< matchContourDistance << "<- HDist:" << dHudist << " Sp:" << distToCentroid << "AreaDist:" << abs(tArea - sArea) << "LengthDist:" << abs(tLength - sLength);

//                   //Reject match 0 in case contour is not actually there
//                   //if (matchContourDistance < gi_ThresholdMatching)
//                        bContourfound = true;
//               }
//           }
//       }


//   if (!bContourfound)
//   {
//       std::cerr << "Failed,Closest Contour :" << idxContour << " d:" << mindistToCentroid << std::endl;
//       idxContour = -1;
//   }
//      //qDebug() << "-------Got best " <<  idxContour << " D:"<< mindistToCentroid;

//   assert(idxContour < (int)contours.size());

//   return idxContour;
//}

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
/// \brief detectZfishFeatures - Used to create geometric representations of main zebrafish Features : Eyes, Body, tail
/// these are saved as point arrays on which angles and other measurements can be obtained
/// \param maskedGrayImg - IMage Masked so only fish is being shown Showing
/// \param imgFishHeadSeg initialized canvas for head
/// \return
///
/// // \todo Optimize by re using fish contours already obtained in enhance fish mask
void detectZfishFeatures(MainWindow& window_main,const cv::Mat& fullImgIn,cv::Mat& fullImgOut,cv::Mat& imgFishHeadSeg,cv::Mat& outimgFishHeadProcessed, cv::Mat& maskedfishImg_gray, std::vector<std::vector<cv::Point> >& contours_body,std::vector<cv::Vec4i>& hierarchy_body)
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
            if (fish == 0 ) //|| fish->inactiveFrames > 1
                continue;

          //Draw A general Region Where the FIsh Is located,
          cv::Point centre = fish->ptRotCentre; //top_left + rotCentre;
          //cv::Point centroid = fish->ptRotCentre ; // cv::Point2f(fish->track->centroid.x,fish->track->centroid.y);
          cv::Point pBound1 = cv::Point(max(0,min(maskedImg_gray.cols,centre.x-gFishBoundBoxSize)), max(0,min(maskedImg_gray.rows,centre.y-gFishBoundBoxSize)));
          cv::Point pBound2 = cv::Point(max(0,min(maskedImg_gray.cols,centre.x+gFishBoundBoxSize)), max(0,min(maskedImg_gray.rows,centre.y+gFishBoundBoxSize)));

          cv::Rect rectFish(pBound1,pBound2);

          //cv::rectangle(fullImgOut,rectFish,CV_RGB(20,200,150),2); //Identify Fish Region Bound In Cyan Square
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
          /// Ellipsoid Vectors detected for each eye
           tEllipsoids vellLeft;
           tEllipsoids vellRight;



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
              imgFishHead           = imgFishAnterior_Norm(rectFishHeadBound);

              /// EYE DETECTION Report Results to Output Frame //
              /// Returns imgFishHeadProcessed Upsampled with ellipses drawns, and imgFishHeadSeg - the processed edges img used
              /// to detect the eyes
              ///
              int ret = 0;
              std::stringstream ss;
              ret = detectEyeEllipses(imgFishHead,vellLeft,vellRight,imgFishHeadSeg,imgFishHeadProcessed);

              if ((ret < 2 | gUserReward < 0) )
              {
                fish->nFailedEyeDetectionCount++;
                if ((fish->nFailedEyeDetectionCount % 100)==0)
                {
                    ss << "-#" << fish->nFailedEyeDetectionCount << " Eye Detection Error - Check Threshold";
                    window_main.LogEvent(QString::fromStdString(ss.str()));
                }
              }

              //Copy Detected Ellipse Frame To The Output Frame
              if (imgFishHeadProcessed.u)
                  imgFishHeadProcessed.copyTo(outimgFishHeadProcessed);

              // Debug test //
              //if (gthresEyeSeg < 0)
              //    gUserReward = -500;

              /// Pass detected Ellipses to Update the fish model's Eye State //
              /// Make Fit score count very little
              double fitScoreReward = 0.0000*fish->updateEyeState(vellLeft,vellRight)+ gUserReward;
              gUserReward = 0.0; //Reset User Provided Rewards
              //qDebug() << "R:" << fitScoreReward;
              tEyeDetectorState current_eyeState = pRLEye->getCurrentState();
              current_eyeState.iSegThres1        = gthresEyeSeg; //Update to what the environment state is
              current_eyeState.iDSegThres2       = gthresEyeSeg-gthresEyeSegL;
              current_eyeState.setVergenceState( fish->leftEyeTheta - fish->rightEyeTheta);
              //Test - Add Target Vector
              //fitScoreReward += -10.0*abs(40.25 - fish->leftEyeTheta) -10.0*abs(34.7 - fish->rightEyeTheta);
              pRLEye->UpdateStateValue(current_eyeState,fitScoreReward); //Tell RL that we moved state so it calc value and updates internal state
              pRLEye->setCurrentState(current_eyeState);
              /// \note RL Is disabled - No State changes
              //tEyeDetectorState new_eyeState = pRLEye->DrawNextAction(current_eyeState); //RL choose next Action /
              //gthresEyeSeg = new_eyeState.iSegThres1; //Action Sets a partial state -> Update threshold as indicated by algorithm
              //gthresEyeSegL = new_eyeState.iSegThres1-new_eyeState.iDSegThres2;
              pwindow_main->UpdateSpinBoxToValue();

              ///  Print Out Values //
              /// \todo Figure out Why/how is it that nan Values Appeared in Output File : NA Values in ./Tracked07-12-17/LiveFed/Empty//AutoSet420fps_07-12-17_WTLiveFed4Empty_286_005_tracks_2.csv
              /// \todo Move this to specialized Function Like @renderFrameText
              ///If Both Eyes Detected Then Print Vergence Angle
              ss.precision(3);
              cv::Scalar colTxt;
              if (fitScoreReward > 100)
                  colTxt = CV_RGB(250,250,50);
              else
                 colTxt = CV_RGB(250,20,20);
              {
                  ss.str(""); //Empty String
                  ss << "L:" << fish->leftEyeTheta;
                  cv::putText(fullImgOut,ss.str(),cv::Point(rect_pasteregion.br().x-45,rect_pasteregion.br().y+10),CV_FONT_NORMAL,0.4,CV_RGB(250,250,0),1 );
                  ss.str(""); //Empty String
                  ss << "R:"  << fish->rightEyeTheta;
                  cv::putText(fullImgOut,ss.str(),cv::Point(rect_pasteregion.br().x-45,rect_pasteregion.br().y+25),CV_FONT_NORMAL,0.4,CV_RGB(250,250,0),1 );
                  ss.str(""); //Empty String
                  ss << "V:"  << ((int)((fish->leftEyeTheta - fish->rightEyeTheta)*10)) /10.0; //Store to global variable
                  cv::putText(fullImgOut,ss.str(),cv::Point(rect_pasteregion.br().x-45,rect_pasteregion.br().y+40),CV_FONT_NORMAL,0.4,CV_RGB(250,250,0),1 );
              }


              //ss.str(""); //Empty String

              //Check If Too Many Eye Detection Failures - Then Switch Template
              if (fish->nFailedEyeDetectionCount > 40)
              {
                    //fish->idxTemplateRow = iLastKnownGoodTemplateRow = (rand() % static_cast<int>(gnumberOfTemplatesInCache - 0 + 1));//Start From RANDOM rOW On Next Search
                    //pwindow_main->LogEvent(QString("[warning] Too Many Eye detection Failures - Change Template Randomly to :" + QString::number(iLastKnownGoodTemplateRow)));
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
                fish->fitSpineToContour2(maskedImg_gray,contours_body,0,idxFish);
#endif

               /// Use Contour Variational Fitting to distance from spine - Adjusts spine segment length to tail contour length
               /// \note If done on all frames is converges on Local Minima where the tail is fit in the body contour.
#ifdef _USEPERIODICSPINETOCONTOUR_TEST
               if (pwindow_main->nFrame == uiStartFrame || (pwindow_main->nFrame % 4  ) == 0)//((uint)gfVidfps/4)
               {
                   int idxFish = findMatchingContour(contours_body,hierarchy_body,centre,2);
                   if (idxFish>=0)
                   {
                        double err_sp0 = fish->fitSpineToContour2(maskedImg_gray,contours_body,0,idxFish);
                   }
                   gFishTailSpineSegmentLength <- fish->c_spineSegL;
                   pwindow_main->UpdateTailSegSizeSpinBox(fish->c_spineSegL);
                   //qDebug() << "Spine Tail Fit Error :" << fish->lastTailFitError;
               }
               /// Main Method Uses Pixel Intensity //
               fish->fitSpineToIntensity(maskedfishFeature_blur,gFitTailIntensityScanAngleDeg);
               fish->drawSpine(fullImgOut);

               //If Convergece TimedOut Then likely the fit is stuck with High Residual and no gradient
               //Best To reset Spine and Start Over Next Time
               /// \todo Make this parameter threshold formal
               if (abs(fish->lastTailFitError) > fish->c_fitErrorPerContourPoint)
               {
                   gFishTailSpineSegmentLength = gc_FishTailSpineSegmentLength_init;
                   fish->c_spineSegL = gFishTailSpineSegmentLength ;
                   pwindow_main->UpdateTailSegSizeSpinBox(fish->c_spineSegL);

                   //fish->resetSpine(); //No Solution Found So Reset
                   //pwindow_main->LogEvent("[info] Reset Spine");

                   pwindow_main->LogEvent(QString("[warning] lastTailFitError ") + QString::number(fish->lastTailFitError) + QString(" > c_fitErrorPerContourPoint") );
//                   int idxFish = findMatchingContour(contours_body,hierarchy_body,centre,2);
//                   double err_sp0 = fish->fitSpineToContour2(maskedImg_gray,contours_body,0,idxFish);
//                   pwindow_main->LogEvent(QString("[info] new lastTailFitError ") + QString::number(fish->lastTailFitError) + QString(" > c_fitErrorPerContourPoint") );
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




