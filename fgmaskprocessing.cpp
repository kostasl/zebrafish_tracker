///*
/// Kostas Lagogiannis 2018
/// File Contains Auxilary functions used to define and clarify Fish And Food Masks
///
///
///*

// #include <opencv2\opencv.hpp>

#include <GUI/mainwindow.h>
#include <fgmaskprocessing.h>
#include <larvatrack.h>

///Curve Smoothing and Matching
#include <CSS/CurveCSS.h>

#include <opencv2/core/ocl.hpp> //For setting setUseOpenCL
/// CUDA //
/// #include <opencv2/opencv_modules.hpp> //THe Cuda Defines are in here
#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    #include <opencv2/core/cuda.hpp>
    #include "opencv2/cudaimgproc.hpp"
    #include "opencv2/cudaarithm.hpp"
    #include <opencv2/photo/cuda.hpp>
    #include <opencv2/core/cuda_types.hpp>
#endif


//const int gcFishContourSize         = ZTF_FISHCONTOURSIZE;


//extern int g_Segthresh;
extern cv::Mat kernelOpen;
extern cv::Mat kernelClose;
extern cv::Mat kernelDilateMOGMask;
extern cv::Mat kernelOpenfish;

extern cv::Point gptHead; ///\todo remove this global var hack

extern double dBGMaskAccumulateSpeed;
// Gaussian Curve Smoothing Kernels For fish Contour//
extern std::vector<double> gGaussian,dgGaussian,d2gGaussian; //These Are init. in main


extern trackerState gTrackerState;
//extern const double dLearningRate; //Learning Rate During Initial BG Modelling done over MOGhistory frames
//extern const double dLearningRateNominal;
//extern double gdMOGBGRatio;
//When Running Multiple Threads That Use BG Substractor - An SEGFault is hit in OpenCL
extern cv::Ptr<cv::BackgroundSubtractorMOG2> pMOG2; //MOG2 Background subtractor

//extern bool bRemovePixelNoise;
//extern bool bUseBGModelling;
//extern bool bUseGPU;
//extern bool bUseOpenCL;
//extern bool bshowMask;

//extern ltROIlist vRoi;
//extern cv::Point ptROI1;
//extern cv::Point ptROI2; //This Default Value Is later Modified
//extern cv::Size gszTemplateImg; //Used For ROI size

extern MainWindow* pwindow_main;

extern fishModels vfishmodels;



#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    extern cv::cuda::GpuMat dframe_gray;
    extern cv::cuda::GpuMat dframe_mask; //Passed to MOG Cuda
    extern cv::cuda::GpuMat dframe_thres; // Used In Mask Enhancement
    extern cv::Ptr<cv::cuda::TemplateMatching> gpu_MatchAlg;// For Template Matching
    extern Ptr<cuda::Filter> gpu_DilateFilter;

#endif

/*// Example Of Mean Image
Mat3b getMean(const vector<Mat3b>& images)
{
    if (images.empty()) return Mat3b();

    // Create a 0 initialized image to use as accumulator
    cv::Mat m(images[0].rows, images[0].cols, CV_64FC3);
    m.setTo(Scalar(0,0,0,0));

    // Use a temp image to hold the conversion of each input image to CV_64FC3
    // This will be allocated just the first time, since all your images have
    // the same size.
    Mat temp;
    for (int i = 0; i < images.size(); ++i)
    {
        // Convert the input images to CV_64FC3 ...
        images[i].convertTo(temp, CV_64FC3);

        // ... so you can accumulate
        m += temp;
    }

    // Convert back to CV_8UC3 type, applying the division to get the actual mean
    m.convertTo(m, CV_8U, 1. / images.size());
    return m;
}
*/

extern QElapsedTimer gTimer;
extern bool bExiting;//Exit Flag



///*
///Create FG Model Image - Since target objects can be/will be are moving from the 1st frame, we need a statistical model
/// of the BG precalculated
//// Uses A Mean IMage Approach to detect Stationary Objects that should be ignored when tracking
/// //Make Sure fgMask Is empty On 1st Call
///
unsigned int getBGModelFromVideo(cv::Mat& bgMask,MainWindow& window_main,QString videoFilename,QString outFileCSV,unsigned int MOGhistoryLength)
{

        cv::Mat frame,frame_gray;

        const int startFrameCount   = 1; //Start Modelling From THe Start
        unsigned int nFrame         = 1; //Current Frame Number
        char keyboard               = 0;

        cv::Mat bgAcc;
        //std::clog << gTimer.elapsed()/60000.0 << " Starting Background Model processing..." << std::endl;
        window_main.LogEvent(" Starting Stat Pixel from "+ QString::number(MOGhistoryLength)+ " frames Background Model processing:" + videoFilename);
        //create the capture object
        cv::VideoCapture capture(videoFilename.toStdString());

        if(!capture.isOpened())
        {
            //error in opening the video input
            std::cerr <<  gTimer.elapsed()/60000.0 << " Unable to open video file: " << videoFilename.toStdString() << std::endl;
            std::exit(EXIT_FAILURE);
        }

        //Get Length of Video
        uint totFrames = capture.get(cv::CAP_PROP_FRAME_COUNT);
        //Do Not Overrun Vid Length In BG COmputation
        uint uiStopFrame = totFrames ; //std::min(totFrames,MOGhistoryLength);
        //read input data. ESC or 'q' for quitting
        uint uiLearnedFrames = 0;
        uint skipFrames = uiStopFrame/ MOGhistoryLength;
        window_main.setTotalFrames(uiStopFrame); //Set To BG Processing REgion

        //  Check If it contains no Frames And Exit
        if (uiStopFrame < 2)
        {
            pwindow_main->LogEvent("[ERROR] This Video File is empty ");
            capture.release();
            return 0;
        }

        while( !gTrackerState.bExiting && (char)keyboard != 27 && nFrame < (uint) uiStopFrame && uiLearnedFrames < MOGhistoryLength)
        {
            uiLearnedFrames++;
            //read the current frame
            if(!capture.read(frame))
            {
                if (nFrame == startFrameCount)
                {
                    std::cerr << gTimer.elapsed()/60000.0 <<  " Unable to read first frame." << std::endl;
                    nFrame = 0; //Signals To caller that video could not be loaded.
                    exit(EXIT_FAILURE);
                }
                else
                {
                   std::cerr << gTimer.elapsed()/60000.0 << ". Unable to read next frame. So this video Is done <<<<<<<" << std::endl;
                   std::clog << gTimer.elapsed()/60000.0 << ". " << nFrame << " frames of Video processed. Move on to next " <<std::endl;
                  //  break;
                   continue;
               }
            } //If Failed to Grab frame
            else
            {//Frame Grabbed - Process It

                assert(!frame.empty());
                //Get Frame Position From Vid Sam
                nFrame = capture.get(cv::CAP_PROP_POS_FRAMES) + startFrameCount;
                window_main.nFrame = nFrame; //Update Window
                window_main.tickProgress();
//              Check If FG Mask Has Been Created - And Make A new One
               if (bgMask.cols == 0)
               {
                    bgMask = cv::Mat::zeros(frame.rows,frame.cols,CV_8UC1);
                    // Add Roi To Mask Otherwise Make On Based oN ROI
//                        cv::circle(bgMask,newROI.centre,newROI.radius,CV_RGB(255,255,255),-1);
               }
               if (bgAcc.empty()) //Make EMpty Mask
                    bgAcc = cv::Mat::zeros(frame.rows,frame.cols,CV_32FC(bgMask.channels()) ); //AccumWeight Result needs to be CV_32FC

               frame.copyTo(frame,bgMask);
               cv::cvtColor( frame, frame_gray, cv::COLOR_BGR2GRAY);

               updateBGFrame(frame_gray, bgAcc, nFrame, MOGhistoryLength);
            }//Frame Grabbed


           checkPauseRun(&window_main,keyboard,nFrame);

           //Jump To Next Frame To Learn
           capture.set(cv::CAP_PROP_POS_FRAMES, nFrame+ skipFrames); //Move To Next So We Take MOGHistory Samples From the Video

        } //main While loop

        //Remove Low Values
        double uiMaxVal,uiMinVal;
        // Threshold Accumulated Mask For Stationary Objects
        ///Find Max Value,this should belong to stationary objects, and Use it as a relative measure to detect BG Objects
        cv::minMaxLoc(bgAcc,&uiMinVal,&uiMaxVal,0,nullptr);

        bgAcc.convertTo(bgMask,CV_8UC1);
        int thres = cv::threshold(bgMask,bgMask,uiMaxVal*0.05,255,cv::THRESH_BINARY | cv::THRESH_OTSU); //All; Above 33% of Max are Stationary
        pwindow_main->LogEvent("Static Food Mask theshold at " + QString::number(thres));

        if (gTrackerState.bStaticBGMaskRemove & gTrackerState.bshowMask)
           cv::imshow("Accumulated BG Model Thresholded",bgMask);

        cv::morphologyEx(bgMask,bgMask, cv::MORPH_CLOSE, kernelDilateMOGMask,cv::Point(-1,-1),1);

        pMOG2->getBackgroundImage(gframeBGImage);
        pwindow_main->showVideoFrame(gframeBGImage,nFrame);

#if defined(DEBUG)
        cv::imshow("BGImage",gframeBGImage);
#endif
        //delete capture object
        capture.release();

        // Save The background image



        //std::clog << gTimer.elapsed()/60000.0 << " Background Processing  loop. Finished" << std::endl;
        window_main.LogEvent(" Background Processing  loop. Finished");


      return nFrame;
} ///trackImageSequencefile




///
/// \brief processMasks Updates BG Model, Can Filter Pixel noise from frame_gray, Adds threshold and BGModel Masks, and it will also mask the BGModel-Mask with the Static mask if bStaticAccumulatedBGMaskRemove  Option Flag is set,
/// if bRemovePixelNoise, then Noise filtering is used to improve FG segmentation (causes slowdown)
/// \returns bgMaskInOut - the combined FG object masks depending on user options
/// \param frame_gray //Current greyScale Frame - Noise May be Removed If filtering Is Set To On
/// \param bgStaticMaskInOut The mask provided to processFrame, Includes Static Objects and ROI Region - if bStaticAccumulatedBGMaskRemove is T, then static mask is combined
/// \param fgFrameOut Returns the Image Of the FG objects only - Can include Fish and moving Prey
/// \param dLearningRate - The speed which MOG learns FG Objects -Set to slow when tracking so fish does not fade when stationary (will be picked up from Threshold mask)
///
void extractFGMask(cv::Mat& frameImg_gray,cv::Mat fgStaticMaskIn,cv::Mat& fgMaskInOut,cv::Mat& fgFrameOut,double dLearningRate)
{
 const int max_thresh = 255;

  cv::Mat threshold_output;

 ////  Obtain FG IMage: Substract MOG Extracted BG Image //
 if (gTrackerState.bUseBGModelling && !gframeBGImage.empty()) {
    fgFrameOut = frameImg_gray - gframeBGImage; //Remove BG Image
 }else{
    std::clog << " [Warning] Missing background model image " << std::endl;
    fgFrameOut  = frameImg_gray;
 }

 //If We are during the Static Mask Accumulation phase ie (fgStaticMaskIn.type() !=  CV_8U) then Produce the Threshold Image
 // The threshold mask is blended with OR with the MOG mask.
 cv::threshold( fgFrameOut, threshold_output, gTrackerState.g_FGSegthresh, max_thresh, cv::THRESH_BINARY ); // Log Threshold Image + cv::THRESH_OTSU

/// \todo: The CUDA code below Needs to Be Revised so it matches the Non-Cuda Process
#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
         if (bUseGPU)
         {
             dframe_gray.upload(frame_gray);

             if (gTrackerState.bRemovePixelNoise)
             {        ///Remove Pixel Noise
                 ///* src – Input 8-bit 1-channel, 2-channel or 3-channel image.         ///        dst – Output image with the same size and type as src .         ///       templateWindowSize – Size in pixels of the template patch that is used to compute weights. Should be odd. Recommended value 7 pixels         ////        searchWindowSize – Size in pixels of the window that is used to compute weighted average for given pixel. Should be odd. Affect performance linearly: greater searchWindowsSize - greater denoising time. Recommended value 21 pixels         ///        h – Parameter regulating filter strength. Big h value perfectly removes noise but also removes image details, smaller h value preserves details but also preserves some noise         ///
                cv::cuda::fastNlMeansDenoising(dframe_gray, dframe_gray,2.0, 21,7);
                dframe_gray.download(frame_gray);
             }
             if (bUseBGModelling)
             {
                 try{
                        pMOG2->apply(dframe_gray,dframe_mask,gTrackerState.dLearningRate);
                    dframe_mask.download(gTrackerState.mMOGMask);
                 }catch(...)
                 {
                     pwindow_main->LogEvent("[Error] CUDA MOG2 failed");
                     //cv::ocl::setUseOpenCL(false); //When Running Multiple Threads That Use BG Substractor - An SEGFault is hit in OpenCL
                 }
             } //BGMOdel
           }//Use GPU
         else{ //Note This Is the Same Code as In THe NO USE_CUDA Case
             if (gTrackerState.bRemovePixelNoise)
                     cv::fastNlMeansDenoising(frame_gray, frame_gray,2.0,7, 21);
                   //Check If BG Ratio Changed
             if (gTrackerState.bUseBGModelling)
             {
               try{
                   pMOG2->apply(frame_gray,gTrackerState.mMOGMask,gTrackerState.dLearningRate);
               }catch(...)
               {
                   std::clog << "MOG2 apply failed, probably multiple threads using OCL, switching OFF" << std::endl;
                   pwindow_main->LogEvent("[Error] MOG2 failed, probably multiple threads using OCL, switching OFF");
                   cv::ocl::setUseOpenCL(false); //When Running Multiple Threads That Use BG Substractor - An SEGFault is hit in OpenCL
                   bUseOpenCL = false;
               }
             }//BGModel
         }
#else //// NON GPU VERSION
      // Remove Pixel Noise
      if (gTrackerState.bRemovePixelNoise)
              cv::fastNlMeansDenoising(frameImg_gray, frameImg_gray,2.0,7, 21);
      //Update BG Model - If not stuck on same frame by being paused
      if (gTrackerState.bUseBGModelling)
      {
        try{
            pMOG2->apply(frameImg_gray,gTrackerState.mMOGMask,dLearningRate);
            //
        }catch(...)
        {
        //##With OpenCL Support in OPENCV a Runtime Assertion Error Can occur /
        //In That case make OpenCV With No CUDA or OPENCL support
        //Ex: cmake -D CMAKE_BUILD_TYPE=RELEASE -D WITH_CUDA=OFF  -D WITH_OPENCL=OFF -D WITH_OPENCLAMDFFT=OFF -D WITH_OPENCLAMDBLAS=OFF -D CMAKE_INSTALL_PREFIX=/usr/local
        //A runtime Work Around Is given Here:
            std::clog << "MOG2 apply failed, probably multiple threads using OCL, switching OFF" << std::endl;
            pwindow_main->LogEvent("[Error] MOG2 failed, probably multiple threads using OCL, switching OFF");
            cv::ocl::setUseOpenCL(false); //When Running Multiple Threads That Use BG Substractor - An SEGFault is hit in OpenCL
        }
      }

#endif

      if (gTrackerState.bUseBGModelling ) //MOG Mask Exists
      {
            //No Static Mask - But Combine With threshold So MOG Ghosts Are Erased
           ///TODO make into and mask , and isolate threshold mask around fish
           cv::bitwise_and(gTrackerState.mMOGMask,threshold_output,fgMaskInOut); // Combine with Threshold So We have some control of output

           //Use MOG mask Not Thresholded OUtput - Improves Boundary Removal
           //fgMOGMask.copyTo(fgMaskInOut);

            //Combine Static Mask and Remove Stationary Learned Pixels From Mask If Option Is Set
            if (gTrackerState.bStaticBGMaskRemove && !fgStaticMaskIn.empty() && fgStaticMaskIn.type() == CV_8U)//Although bgMask Init To zero, it may appear empty here!
                  cv::bitwise_and(fgStaticMaskIn,fgMaskInOut,fgMaskInOut); //Only On Non Stationary pixels - Ie Fish Dissapears At boundary
       }else{  // No MOG Mask Exists so simply Use thresh Only to Detect FG Fish Mask Detect
          //Returning The thresholded image is only required when No BGMask Exists
         if (gTrackerState.bStaticBGMaskRemove && !fgStaticMaskIn.empty() && fgStaticMaskIn.type() == CV_8U)
            cv::bitwise_and(fgStaticMaskIn,threshold_output,fgMaskInOut);
         else //No Static Mask here and no BG Model, so simply return the thresholded one
             threshold_output.copyTo(fgMaskInOut);
      } //Threshold Only Available

      // REMOVE Items Outside ROI


      // Show Masks for Debuging purposes
      if (gTrackerState.bshowMask && !fgStaticMaskIn.empty())
        cv::imshow("StaticMask",fgStaticMaskIn);
      if (gTrackerState.bshowMask && !gTrackerState.mMOGMask.empty())
        cv::imshow("MOGMask",gTrackerState.mMOGMask);
     if (gTrackerState.bshowMask && !threshold_output.empty())
        cv::imshow("Threshold Out",threshold_output);

} //END PROCESSMASKS

/// \brief Find Contours Max Inflection point/Sharpest point
/// Traces the Point/Idx of the sharpest contour point, then it Smooths and Simplifies the curve (contour)
/// Adding the Sharp point back into the curve
/// \returns idx Of Sharp point in the curve provided
int getMaxInflectionAndSmoothedContour(const cv::Mat& frameImg, cv::Mat& fgMask,std::vector<cv::Point>& curve)
{

    ///// SMOOTH COntours /////
    vector<double> curvex,curvey,smoothx,smoothy,resampledcurveX,resampledcurveY ;
    PolyLineSplit(curve,curvex,curvey);

    std::vector<double> X,XX,Y,YY;
    std::vector<double> dXY;

    getdXcurve(curvex,gTrackerState.dGaussContourKernelSigma ,smoothx,X,XX,gTrackerState.gGaussian,gTrackerState.dgGaussian,gTrackerState.d2gGaussian,false);
    getdXcurve(curvey,gTrackerState.dGaussContourKernelSigma,smoothy,Y,YY,gTrackerState.gGaussian,gTrackerState.dgGaussian,gTrackerState.d2gGaussian,false);

/// Curve Signature Analysis ///
/*    //Finds Inwards Curvature points that exceed a threshold, defined as the maximal curvature points
    vector<int> vidxMax = ComputeCSSImageMaximas(curvex,curvey,smoothx,smoothy);
    vector<vector<Point> > contours(1);
    PolyLineMerge(contours[0], smoothx, smoothy);
    cv::Mat contourimg;
    fgMask.copyTo(contourimg);
    cv::drawContours(contourimg, contours, 0, Scalar(255,255,255),1, cv::LINE_8);
    for (vector<int>::iterator itr = vidxMax.begin(); itr!=vidxMax.end(); ++itr) {
        cv::circle(contourimg, contours[0][*itr],4,CV_RGB(255,255,255),cv::FILLED);
    }
    cv::imshow("contour",contourimg);
*/
/// End OF Curve Analysis //

    dXY.resize(X.size());

    /// Find Tail As POint Of Maximum Curvature dXY
    int idxMin2 = 0;
    int idxMin  = 0;
    double maxVal=0.0;
    double minVal=10000.0;
    cv::Point ptSharp,ptSharp2,ptHead,ptHead2,ptTail;

    for (int j=0; j<(int)X.size(); j++) {
       dXY[j] = (X[j]*X[j] + Y[j]*Y[j]);
       maxVal = dXY[j];

       if (dXY[j] < minVal) //Detect Tail
       {
           idxMin2 = idxMin;//Save As 2nd Smallest
           idxMin = j;
           minVal = dXY[j];
       }
    }//for loop

    ptSharp = curve[idxMin]; //Most Likely Tail Point
    ptSharp2 = curve[idxMin2]; //This Could Be Head

    //Only Simplify Points if contour Curve Too long
    if ((int)curve.size() > (2*gTrackerState.gcFishContourSize))
    {
        ResampleCurve(smoothx,smoothy,resampledcurveX,resampledcurveY,gTrackerState.gcFishContourSize,false);
        PolyLineMerge(curve,resampledcurveX,resampledcurveY); //Remerge
    }


    //Find Where Tail Point Is In the Resampled (Reduced) Contour
    int idxTail = findIndexClosesttoPoint(curve,ptSharp);

    if (idxTail >= 0)
    {
        //Put Tail Back to Curve In CAse it Got Smoothed Out
        std::vector<cv::Point>::iterator it = curve.begin();
        it += idxTail;
        curve.insert(it,ptSharp);
    }

    return(idxTail);
}
/// END OF Max INflection Tail Detect ///

/// \brief Returns the index of the point furthest away from the provided tail point idx on the curve
/// That furthest from the tail position is likely the larva's head
int findAntipodePointinContour(int idxTail, std::vector<cv::Point>& curve,cv::Point ptCentroid, cv::Point2f& ptHead,cv::Point2f& ptTail)
{
    cv::Point ptHead2;
    int retIdx;
    assert(curve.size() > 0);

    // Copy with a reshuffling of TailPoint to the zero index
    vector<cv::Point> vcurveR( curve.begin() + idxTail ,curve.end());
    vcurveR.insert(vcurveR.end(),curve.begin(),curve.begin() + std::max(0,idxTail) );

    /// Find Head-Tail Point As Point Of Maximum Arc Length //
    //1st Approach is the noddy Way, n!
    int maxLen = 0;
    int lenA,lenB; //Arc Lenght ClockWise, And AntiClockWise
    int idxA,idxB,idxT;
    idxA = idxB = idxT = lenA = lenB = maxLen = 0;

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
    /// CHECK Here For BUG
    assert( idxTail < curve.size() );
    assert( idxB < vcurveR.size() );
    //Check Which curve point is closest to the TailInflection Candidate - Initially Mark As Tail
    if (norm(curve[idxTail]-vcurveR[idxA]) < norm(curve[idxTail]-vcurveR[idxB]) )
    {
        ptTail = vcurveR[idxA];
        ptHead = vcurveR[idxB];
        ptHead2 = vcurveR[idxB-1];
        retIdx = idxB;
    }
    else
    {
        ptTail = vcurveR[idxB];
        ptHead = vcurveR[idxA];
        ptHead2 = vcurveR[idxA-1];
        retIdx = idxA;
    }
     //Verify Head Point is closest to the Centroid (COM) than tail / Otherwise Switch
    if (norm((cv::Point)ptTail-ptCentroid)  < norm((cv::Point)ptHead-ptCentroid))
    {
        cv::Point temp = ptTail;
        ptTail = ptHead;
        ptHead = temp;
    }

    //Replace Old Curve
    curve = vcurveR;

    return retIdx; //Return Head Idx in Curve
}

/// \brief handles the processing of Prey Item Mask such that prey detection can be improved.
void getPreyMask(const cv::Mat& frameImg, cv::Mat& fgMask,cv::Mat& outFoodMask)
{

    /// Create Thresholded Image Mask from FG image ///
    cv::threshold( frameImg, outFoodMask, gTrackerState.g_SegFoodThesMin, 255, cv::THRESH_BINARY ); // Log Threshold Image + cv::THRESH_OTSU

    /// TODO: Change this As it Gets stuck With Excluding
    if (gTrackerState.bUseBGModelling && !fgMask.empty()) //We Have a (MOG) Model In fgMask - So Remove those Stationary Pixels
        /// \note frameImg may already have FG mask applied to it
        //bitwise_and(outFoodMask,gTrackerState.mMOGMask,outFoodMask); //Apply MOG mask Only
        bitwise_and(outFoodMask,fgMask,outFoodMask);
    ///Removed because DEnse prey conditions Result in Large Blobs
    //Shrink Dilate prey so as to improve tracking and remove Noise
    cv::morphologyEx(outFoodMask,outFoodMask,cv::MORPH_OPEN,kernelDilateMOGMask,cv::Point(-1,-1),1); //
    cv::morphologyEx(outFoodMask,outFoodMask,cv::MORPH_CLOSE,kernelDilateMOGMask,cv::Point(-1,-1),1); //

}

/// \brief Calculate Angle Between Pt 1 and pt 2 in degrees
double ptangle_deg(const Point& v1, const Point& v2)
{
    double cosAngle = v1.dot(v2) / (cv::norm(v1) * cv::norm(v2));
    if (cosAngle > 1.0)
        return 0.0;
    else if (cosAngle < -1.0)
        return CV_PI;
    return std::acos(cosAngle)* 180 / CV_PI;
}


/// \brief handles the processing of Fish Item Mask such that larval detection can be improved.
///  uses optic flow info to detect if blob belongs to a fishmodel so as to classify it correctly
std::vector<std::vector<cv::Point> > getFishMask(const cv::Mat& frameImg_grey,const cv::Mat& outUserFrame, cv::Mat& fgMask, cv::Mat& outFishMask, zftblobs& ptFishblobs)
{
    cv::Mat mask_fnetScore,imgFishAnterior_NetNorm; //For FishNet Detect
    cv::Mat fgEdgeMask;
    std::vector<std::vector<cv::Point> > vFilteredFishbodycontours;
    std::vector<std::vector<cv::Point> > vFilteredFishbodycontours_classified;
    zftblobs vFishKeypoints_next;
    cv::Point2f ptHead,ptHead2,ptTail,ptTail2; //Traced Points for Tail and Head
    cv::Point ptSearch; //Helper pt variable
    //Make Empty Canvas to Draw Fish Mask
    outFishMask = cv::Mat::zeros(frameImg_grey.rows,frameImg_grey.cols,CV_8UC1);

    //Shring-Grow -Erase Thin Border Lines
    cv::morphologyEx(fgMask,fgMask,cv::MORPH_DILATE, kernelDilateMOGMask,cv::Point(-1,-1),2);
    //cv::morphologyEx(fgMask,fgMask, cv::MORPH_OPEN, kernelOpenfish,cv::Point(-1,-1),2);
    cv::morphologyEx(fgMask,fgMask, cv::MORPH_CLOSE, kernelOpenfish,cv::Point(-1,-1),3);

    //Make Hollow Mask Directly - Broad Approximate -> Grows outer boundary
    cv::morphologyEx(fgMask,fgEdgeMask, cv::MORPH_GRADIENT, kernelOpenfish,cv::Point(-1,-1),1);

    /// Find contours main Internal and External contour using on Masked Image Showing Fish Outline
    /// //Used RETR_CCOMP that only considers 1 level children hierachy - I use the 1st child to obtain the body contour of the fish
    std::vector<std::vector<cv::Point> > fishbodycontours;
    std::vector<cv::Vec4i> fishbodyhierarchy;

    //Then Use ThresholdImage TO Trace More detailed Contours
    //cv::dilate(threshold_output_COMB,threshold_output_COMB_fish,kernelOpenfish,cv::Point(-1,-1),4);
    assert(!fgEdgeMask.empty());


    ptFishblobs.clear();

    // UnMask Regions Around Existing Fish
    fishModels::iterator ft;
    for ( ft  = vfishmodels.begin(); ft!=vfishmodels.end(); ++ft)
    {
        fishModel* pfish = ft->second;
        //circle(outFishMask,pfish->bodyRotBound.center,50,CV_RGB(255,255,255),cv::FILLED);
       // drawContours(fgEdgeMask,pfish->contour,0,CV_RGB(255,255,255),cv::FILLED);
    }
    if (gTrackerState.bshowMask)
        cv::imshow("FishMask Edge",fgEdgeMask);
    cv::findContours( fgEdgeMask, fishbodycontours,fishbodyhierarchy, cv::RETR_CCOMP,
                      cv::CHAIN_APPROX_SIMPLE , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE

    /* Obtain Optic Flow of fish positions from previous frame - Score Contours as fish based on whether they contain a fish keypoint that moved into the vicinity //
     processFishOpticFlow(fgEdgeMask,gframeLast, vfishmodels, vFishKeypoints_next);
     /// Check if Blob belongs to moving fish - Draw Reveal Mask
        bool  = false;
        for (int k=0; k < vFishKeypoints_next.size();k++)
        {
            if (boundEllipse.boundingRect().contains(vFishKeypoints_next[k].pt ))
            {
                circle(outFishMask,vFishKeypoints_next[k].pt,20,CV_RGB(255,255,255),cv::FILLED);
                bFishBlobFlowed = true;
            }
        }


        if (boundEllipse.size.width/boundEllipse.size.height < 2.5 &&
                boundEllipse.size.height/boundEllipse.size.width < 2.5 &&
                !bFishBlobFlowed)
        {
            if (bFishBlobFlowed)
                qDebug() << "Bounded Shape/Flow Filt.";

            continue;
        }*/


    ///Draw Only the largest contours that should belong to fish
    /// \todo Other Match Shapes Could be used here
    /// \todo Use WaterShed - Let MOG mask Be FG label and then watershed
    //int idxFishContour = -1;
    std::vector<cv::Point> curve; // THe Fish Contour to use for new Mask
    int iHitCount = 0;
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

        // Lets try simple area filter - Assume no large object need to be BG substracted
        uint area  = cv::contourArea(fishbodycontours[kk]);

        cv::Point centroid;

        /// DO BASIC Contour Filtering for Fish Like Features //
        // Area check and then  Find the thresholded Fish Contour std::max(dMeanBlobArea*8,(double)thresh_fishblobarea)
        if (area <  gTrackerState.thresh_minfishblobarea || area > gTrackerState.thresh_maxfishblobarea) //If Contour Is large Enough then Must be fish
            continue; //skip to next Contour

        // Check Centroid is in ROI - Find Barycenter - Centroid
        cv::Moments moments =  cv::moments(fishbodycontours[kk]);
        centroid.x = moments.m10/moments.m00;
        centroid.y = moments.m01/moments.m00;
        //ptHead.x = E_x1; ptHead.y=E_y1;
        //ptTail.x = E_x2; ptTail.y=E_y2;
        //If Contained In ROI
        if (!pointIsInROI(centroid,2)) //
            continue;

        curve = fishbodycontours[kk];

        //If Has enough Points Skip Very Small Curves //If curve is empty then  Small Area Contour will be skipped
        if ((int)curve.size() < gTrackerState.gcFishContourSize/2)
            continue;

        /// Fit Ellipse And Obtain Orientation and Focal Points //
        cv::RotatedRect boundEllipse = cv::fitEllipse(curve);
        //cv::ellipse(outFishMask,boundEllipse,CV_RGB(255,255,255),1,cv::LINE_8);
        // Get  Ellipse Major/Minor Axis
        Point2f vertices[4];
        boundEllipse.points(vertices);
        cv::Point2f  ptEllFocus1 = (vertices[0] + vertices[1])/2.0;
        cv::Point2f  ptEllFocus2 = (vertices[2] + vertices[3])/2.0;

        cv::Point2f  ptEllFocus3 =  (vertices[1] + vertices[2])/2.0;
        cv::Point2f  ptEllFocus4 = (vertices[3] + vertices[0])/2.0;
        //Obtain Major axis  Based on Distances - Assign Head - Tail To Major Axis
        if (cv::norm(ptEllFocus1-ptEllFocus2)<cv::norm(ptEllFocus3-ptEllFocus4)){
           gptHead = ptHead2 = ptHead = ptEllFocus3;
           gptTail = ptTail2 = ptTail = ptEllFocus4;
        }else
        {
           gptHead = ptHead2 = ptHead = ptEllFocus1; //Direction ignorant - Only Orientation is picked up
           gptTail = ptTail2 = ptTail = ptEllFocus2;
        }

        /// Make Blob / Keypoint//
        zftblob kp(boundEllipse.center.x,boundEllipse.center.y,area,(int)(boundEllipse.angle+90)%360);

        if (!pointIsInROI(kp.pt,2)) //2nd ROI Check On Bounded Ellipse
            continue;

        /// \todo Here I could Use Shape Similarity Filtering - Through the CurveCSS header
        //Find Tail Point- As the one with the sharpest Angle
        // Smooth Contour and Get likely Index of Tail point in contour, based on curvature sharpness / And
       ptSearch = gptHead;// kp.pt;
       int idxTail,idxHead;
       idxTail = getMaxInflectionAndSmoothedContour(frameImg_grey, fgMask, curve);
       if (idxTail >= 0 )
       {
           idxHead = findAntipodePointinContour(idxTail,curve,centroid,ptHead2,ptTail2);
           // Classify Keypoint for fish  - Find Best Angle if 1st Pass Fails //
           ptSearch  = (ptHead-ptHead2)/2.0f+ptHead2; //Classifier search focused between candidate Ellipse head point And  Contour Head POint
           kp.pt = ptSearch;
           // Fix Orientation against Head Point detected by classifierkp
           // Carefull Cause sometime the fitted ellipse is too long anteriorly making distance wise appear as the tail
       }

       //Make scanned image Bounds at least gFishBoundBoxSize each side
       boundEllipse.size.height = max(gTrackerState.gFishBoundBoxSize, (int)boundEllipse.size.height);
       boundEllipse.size.width  = max(gTrackerState.gFishBoundBoxSize, (int)boundEllipse.size.width);

       vFilteredFishbodycontours.push_back(curve);

       std::vector<float> fRH(4);  //Score Array For Each Pass

        kp.response = 0;
       //Optimization 1 - Small Area Search Check If kp Belogs to Existing Model Fish - Update position To tracked Head Pt
       for (fishModels::iterator it=vfishmodels.begin(); it!=vfishmodels.end(); ++it)
       {
             fishModel* fish = (*it).second;
             // Blob Is in Existing Fish Contour - Must be pointing to existing model
             if (cv::norm(fish->ptRotCentre - kp.pt) <= 20 )//(fish->bodyRotBound.boundingRect().contains(kp.pt))
             {
               zftblob kpFishBlob = kp; //Search Very Near Last Known Fish Position
               kpFishBlob.pt = fish->ptRotCentre;
               kpFishBlob.response = gTrackerState.fishnet.scoreBlobRegion(frameImg_grey, kpFishBlob,boundEllipse, imgFishAnterior_NetNorm,
                                                                mask_fnetScore,20,2,2, QString::number(iHitCount).toStdString(),false);
              if (kpFishBlob.response > gTrackerState.fishnet_classifier_thres){
                  kp = kpFishBlob;
                  break;
              }
             }
       } //Check for Existing FishModel

       /// Low Match Score Could mean Fish is No longer on the same position -
       ///  Scan Region Around Blob to Find precise Head point - using Classifier
       if (kp.response < gTrackerState.fishnet_classifier_thres)
       {
               /// PYRAMID Detection
               //1st Pass - blob position is updated to Detected Position
               fRH[0] = gTrackerState.fishnet.scoreBlobRegion(frameImg_grey, kp,boundEllipse, imgFishAnterior_NetNorm,
                                                                mask_fnetScore,100,10,10, QString::number(iHitCount).toStdString(),false);
               if (gTrackerState.bshowDetectorDebugImg)
                    cv::circle(outUserFrame,kp.pt,10,CV_RGB(5,20,5),1,LINE_AA );

               //cv::circle(outFishMask,kp.pt,4,CV_RGB(155,155,155),2);
               //2nd Pass
               if (fRH[0] >= gTrackerState.fishnet_classifier_thres*0.75){
                    fRH[1] = gTrackerState.fishnet.scoreBlobRegion(frameImg_grey, kp,boundEllipse, imgFishAnterior_NetNorm,
                                                                mask_fnetScore,18,3,3, QString::number(iHitCount).toStdString()+"B",false);
                    if (gTrackerState.bshowDetectorDebugImg)
                        cv::circle(outUserFrame,kp.pt,5,CV_RGB(5,80,5),1,LINE_AA );
               }

               else
                    //Show User where Failed Detection Happened
                    cv::circle(outUserFrame,kp.pt,2,CV_RGB(255,5,5),2,LINE_AA );

               //3rd Pass - Extensive Search in Small Region to match best score
               if (fRH[1] >= gTrackerState.fishnet_classifier_thres*0.85)
               {
                    fRH[2] = gTrackerState.fishnet.scoreBlobRegion(frameImg_grey, kp,boundEllipse, imgFishAnterior_NetNorm,
                                                                mask_fnetScore,4,1,1, QString::number(iHitCount).toStdString()+"C",true);
                    if (gTrackerState.bshowDetectorDebugImg)
                        cv::circle(outUserFrame,kp.pt,2,CV_RGB(5,120,5),1,LINE_AA );
               }

               else
                    //Show User LaseFailed Detection Happened
                    cv::circle(outUserFrame,kp.pt,1,CV_RGB(255,5,5),2,LINE_AA );

//               //4th Pass - Extensive Search in Small Region to match best score
//               if (fRH[2] >= gTrackerState.fishnet_classifier_thres*0.95){
//                    fRH[3] = gTrackerState.fishnet.scoreBlobRegion(frameImg_grey, kp,boundEllipse, imgFishAnterior_NetNorm,
//                                                                mask_fnetScore,3,1,1, QString::number(iHitCount).toStdString()+"D",false);

//                    if (gTrackerState.bshowDetectorDebugImg)
//                        cv::circle(outUserFrame,kp.pt,1,CV_RGB(5,200,5),1,LINE_AA );
//               }

//               else
//                    //Show User LaseFailed Detection Happened
//                    cv::circle(outUserFrame,kp.pt,1,CV_RGB(255,5,5),1,LINE_AA );


           }else{ //Match Score was Above threshold
               //Show Where Fish-Contour Match Happened
               cv::circle(outUserFrame,kp.pt,3,CV_RGB(5,5,255),2,LINE_AA );
               //continue; //Next Contour
          } //Check If Key Point Does not belong existing Fish position / and Scan Area to relocate Fish


       // Head Has been located by Classifier - Now Fix Orientation
       //Correct Angle - 0 is Vertical Up
       if (cv::norm(ptHead2-kp.pt)>cv::norm(ptTail2 - kp.pt)){
          gptHead = kp.pt; //Hack To Get Head position to Classify fish image
          gptTail = ptHead2; //On Curve Point
       }else{
          gptHead = kp.pt; //Hack To Get Head position to Classify fish image
          gptTail = ptTail2; //On Curve Point
       }
       kp.angle = (int)(cv::fastAtan2(gptHead.y-gptTail.y,gptHead.x-gptTail.x)+90)%360;

        QString strfRecScore = QString::number(kp.response,'g',3);
        iHitCount++;
        //qDebug() << kp.angle << " (" << kp.pt.x << "," << kp.pt.y << ")" << "R:" << strfRecScore;

        /// Add TO Filtered KP - IF keypoint is still within roi (moved by classifier) and Passes Classifier threshold
        if (kp.response >= gTrackerState.fishnet_classifier_thres &&
                pointIsInROI(kp.pt,2))
        {

            /// After Classifier Identifies Fish Then Add Template Matching - To correct Orientation
            ptSearch = kp.pt;
            if (gTrackerState.bUseTemplateMatching){
                int bestAngle = kp.angle;
                int iTemplRow = gTrackerState.iLastKnownGoodTemplateRow; //Starting Search Point For Template
                int iTemplCol = kp.angle;
                float maxMatchScore = doTemplateMatchAroundPoint(frameImg_grey, kp.pt, iTemplRow, iTemplCol, bestAngle, ptSearch, outUserFrame);//fishblob->response; //  gTrackerState.gTemplateMatchThreshold*1.1;//
                gTrackerState.iLastKnownGoodTemplateRow = iTemplRow;
                //If Template Match Succeeds then Update Blob To Correct Position And Orientation -
                if (maxMatchScore >= gTrackerState.gTemplateMatchThreshold)  //gTrackerState.fishnet_classifier_thres //|| maxMatchScore > 100.0f
                {
                    // Update to template Matched Angle and positio
                    kp.angle            = bestAngle;
                    kp.pt               = ptSearch; //Does not Work Accuratelly
                    //qDebug() << "+Tmpl:" << maxMatchScore;
                }
                else //Fish May Have turned Too Fast And the Restrained Angle Optimization Failed
                { //Search Again Along All Angles
                    iTemplCol = 0;
                    maxMatchScore = doTemplateMatchAroundPoint(frameImg_grey, kp.pt, iTemplRow, iTemplCol, bestAngle, ptSearch, outUserFrame);//fishblob->response; //  gTrackerState.gTemplateMatchThreshold*1.1;//
                    gTrackerState.iLastKnownGoodTemplateRow = iTemplRow;
                    // Update to template Matched Angle and position
                    kp.angle            = bestAngle;
                    kp.pt               = ptSearch; //Does not Work Accuratelly

                }

            }//If template Matching is used

            ptFishblobs.push_back(kp);
            /// \todo Conditionally add this Contour to output if it matches template.
            vFilteredFishbodycontours_classified.push_back(curve);


            ///  COMBINE - DRAW CONTOURS
            ///\bug drawContours produces freezing sometimes
            //Draw New Smoothed One - the idx should be the last one in the vector
            cv::drawContours( outFishMask, vFilteredFishbodycontours_classified, 0, CV_RGB(255,255,255), cv::FILLED); //
            cv::drawContours( outFishMask, vFilteredFishbodycontours_classified, 0, CV_RGB(255,255,255),2,cv::LINE_AA); //
            //cv::drawContours( outUserFrame, vFilteredFishbodycontours_classified, 0, CV_RGB(5,255,5),1,cv::LINE_AA); //


            /// DEBUG - Show imgs
            if (!imgFishAnterior_NetNorm.empty() && gTrackerState.bshowDetectorDebugImg){
                cv::imshow((QString("FishNet Norm ") + QString::number(iHitCount)).toStdString() ,imgFishAnterior_NetNorm);
            } //If ImgAnt Exists

          if (gTrackerState.bTrackedOneFishOnly && vFilteredFishbodycontours_classified.size() > 0)
              break; //Found A blob Exit Loop //Otherwise I need to Cluster Nearby Blobs so I don get double hits
        } //If DNN classifier Success


        //Add Trailing Expansion to the mask- In Case End bit of tail is not showing (ptTail-ptHead)/30+
        cv::circle(outFishMask, gptTail,5,CV_RGB(255,255,255),cv::FILLED);
        cv::circle(outFishMask, gptHead,8,CV_RGB(255,255,255),cv::FILLED);
        cv::circle(outFishMask,ptSearch,3,CV_RGB(255,255,255),2);
        cv::putText(outFishMask,strfRecScore.toStdString(), ptSearch + cv::Point(10,-10), gTrackerState.trackFnt, gTrackerState.trackFntScale ,  CV_RGB(255,255,250));
        //cv::circle(outFishMask, ptHead - (ptHead-centroid)/2,5,CV_RGB(255,255,255),cv::FILLED);
    } //For Each Fish Contour

    //Check If Any Contours Where Classified As Fish
    if (vFilteredFishbodycontours_classified.size() == 0)
    {
        //qDebug() << "No Fish Blob detected.";
        cv::drawContours( outUserFrame, vFilteredFishbodycontours, 0, CV_RGB(255,5,5),1,cv::LINE_AA); //
    }

    return (vFilteredFishbodycontours);

}

///
/// \brief enhanceMask Attempts to separate Prey and Fish Masks
/// It Looks for fish countours and draws them onto the FG mask so as to isolate the Fish and enhance features
/// * Filters contours for area, to pickup fish like ones
/// * Smooths large contour curves and identies maximum curve point as tail (sharp change in curvature)
/// *The opposite (antipode) point to the tail in curve is identified as the head.
/// \param frameImg - Raw Input camera input in Mat - colour or gray -
/// \param fgMask - (IN) FG Mask Image from processMask
/// \param outFishMask -(OUT) Mask Enhanced for Fish Blob Detection
/// \param outFoodMask -(OUT) Enhanced for Food Blob Detection
/// \todo Cross Check Fish Contour With Model Position
/// - Tracker Picks Up Wrong contour Although Template Matching Finds the fish!
/// Note: Should Use MOG Mask for Blob Detect, But . But thresholded IMg For Countour FInding
void enhanceMasks(const cv::Mat& frameImg_grey, cv::Mat& fgMask,cv::Mat& outFishMask,cv::Mat& outFoodMask,cv::Mat& outUserFrame,
                  std::vector<std::vector<cv::Point> >& outfishbodycontours, zftblobs& ptFishblobs)
{

    cv::Mat frameImg_gray_masked;

    std::vector<std::vector<cv::Point> > vFilteredFishbodycontours;
    cv::Point ptHead,ptHead2,ptTail; //Traced Points for Tail and Head

    frameImg_grey.copyTo(frameImg_gray_masked,fgMask);
    //frameImg_gray = frameImg.clone();//frameImg.clone();

    // Check this Again / What is it doing?
    getPreyMask(frameImg_grey,fgMask,outFoodMask);
    gTrackerState.mfgPreyMask = outFoodMask;
    vFilteredFishbodycontours = getFishMask(frameImg_grey,outUserFrame,fgMask,outFishMask,ptFishblobs);
    gTrackerState.mfgFishMask  = outFishMask; //Save MAsk globally
    frameImg_grey.copyTo(gTrackerState.mfgFishFrame,outFishMask);//Apply Mask and Save FG Image Globally
    //Write The fish contour Mask on Food Mask To erase isolated fish Pixels by Using Smoothed Contour
    // Can invert Fish Mask And Apply on to Food Mask
    // Fish Mask Is Only Outline And So it Does not Erase fish
//    cv::Mat fgNonFish;
//    bitwise_not(outFishMask,fgNonFish);
//    bitwise_and(outFoodMask,fgNonFish,outFoodMask);

    outfishbodycontours = vFilteredFishbodycontours;

    if (gTrackerState.bshowMask)
    {
        #if defined(USE_CUDA)
            if (bUseGPU) dframe_thres.download(threshold_output);
        #endif
       cv::imshow("Food Mask",outFoodMask); //Hollow Blobs For Detecting Food
       cv::imshow("Fish Mask",outFishMask);
       if (!fgMask.empty())
           cv::imshow("BG Model",fgMask);
    }

#if defined(_DEBUG)
    cv::imshow("threshold_output_COMB",threshold_output_COMB);
    cv::imshow("Fish Mask",outFishMask);
#endif


}
// End of Enhance Mask

///
/// \brief updateBGFrame Update BG model for a fixed number of frames / Construct Accumulated Model -
/// \callergraph getBGModelFromVideo
/// \param frame
/// \param fgMask
/// \param nFrame
/// \return returns false when limit of updates is reached
///
bool updateBGFrame(cv::Mat& frameImg_gray, cv::Mat& bgAcc, unsigned int nFrame,uint MOGhistory)
{
    cv::Mat bgMaskThresholded;
    zftblobs ptFishblobs;
    std::vector<std::vector<cv::Point> > fishbodycontours;
    std::vector<cv::Vec4i> fishbodyhierarchy;
    bool ret = true;
    //Speed that stationary objects are removed
   // double dblRatioPxChanged    = 0.0;


    // Detect Food at Lower Thresh //
    cv::Mat bgMask,fgFishMask,fgFoodMask,fgFrameImg;

   // cv::equalizeHist( frame, frame );
    //Update MOG,filter pixel noise and Combine Static Mask
    extractFGMask(frameImg_gray,bgAcc,bgMask,fgFrameImg,gTrackerState.dLearningRate); //Applies MOG if bUseBGModelling is on
    ///Enhance Mask With Fish Shape
    enhanceMasks(fgFrameImg,bgMask, fgFishMask, fgFoodMask,frameImg_gray, fishbodycontours, ptFishblobs);

    pwindow_main->showVideoFrame(bgMask,nFrame);
    //Accumulate things that look like food / so we can isolate the stationary ones
    //cv::threshold( fgFrameImg, bgMaskThresholded, gTrackerState.g_FGStaticMaskSegthresh, 255, cv::THRESH_BINARY ); // Log Threshold Image + cv::THRESH_OTSU
    cv::accumulateWeighted(bgMask,bgAcc,gTrackerState.dBGMaskAccumulateSpeed);
    //cv::imshow("bAcc",bgAcc);

    return ret; //If False then tell calling function to stop updating
}





///
/// \brief findMatchingContour Looks for the inner contour in a 2 level hierarchy that matches the point coords
/// \param contours source array in which to search
/// \param hierarchy
/// \param pt - Position around which we are searching
/// \param level - The required hierarchy level description of the contour being searched for
/// - if -1 (default) do not check hierarchy,
/// - 0 only parent contours
/// - 1 Needs to Have a parent
/// - level == 2 Needs to be top Level Contour
/// \return Index of *child*/Leaf contour closest to point
///
int findMatchingContour(std::vector<std::vector<cv::Point> >& contours,
                              std::vector<cv::Vec4i>& hierarchy,
                              cv::Point pt,
                              int level=-1)
{
    int idxContour           = -1;
    bool bContourfound       = false;
    int mindistToCentroid    = +10000; //Start Far
    int distToCentroid       = +10000;
    int matchContourDistance = 10000;


    assert((level > 0) &&  hierarchy.size() ==contours.size() || level ==-1  );


   /// Render Only Countours that contain fish Blob centroid (Only Fish Countour)
   ///Search Through Contours - Draw contours + hull results

   ///Find Contour with Min Distance in shape and space -attach to closest contour
   //In Not found Search Again By distance tp Full Contour
       //Find Closest Contour
       for( int i = 0; i< (int)contours.size(); i++ )
       {

          //Filter According to desired Level
          if (level == 0) /// Only Process Parent Contours
          {
            if (hierarchy[i][3] != -1) // Need to have no parent
               continue;
            if (hierarchy[i][2] == -1)  // Need to have child
                continue;
            assert(hierarchy[hierarchy[i][2]][3] == i ); // check that the parent of the child is this contour i
          }

          if (level == 1) // Only Process Child Contours
          {
              if (hierarchy[i][3] == -1) // Need to have a parent
                  continue;
//                   //Parent should be root
//                   if (hierarchy[hierarchy[i][3]][3] != -1)
//                       continue;
          }

          if (level == 2) // Needs to be top Level Contour
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

       } //Loop through Contours


   if (!bContourfound)
   {
       //std::cerr << "Failed,Closest Contour :" << idxContour << " d:" << mindistToCentroid << std::endl;
       idxContour = -1;
   }
      //qDebug() << "-------Got best " <<  idxContour << " D:"<< mindistToCentroid;

   assert(idxContour < (int)contours.size());

   return idxContour;
}


// Attemts to Get Orientation and centre of an Elonged - FishLike Blob
int getFishBlobCentreAndOrientation(cv::Mat imgFishAnterior,cv::Point2f ptCentre,int Angle,cv::Point2f& ptRevised,int& RevisedAngle)
{
    cv::Mat imgFishAnterior_blur,imgFishAnterior_thres;
    // Correct Orientation and Centre
    std::vector<std::vector<cv::Point> > contours;
    std::vector<cv::Vec4i> hierarchy;
    int iAngleOffset = Angle;
    cv::Point2f ptCentreCorrection = ptCentre;

    cv::GaussianBlur(imgFishAnterior,imgFishAnterior_blur,cv::Size(7,7),5,5);
    cv::adaptiveThreshold(imgFishAnterior_blur, imgFishAnterior_thres,255,cv::ADAPTIVE_THRESH_MEAN_C,cv::THRESH_BINARY,31,0); //Last Param Is const substracted from mean

    /// Find contours //
    cv::findContours( imgFishAnterior_thres, contours, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );


    int idxContour = 0;
    //Should be one contour
    if (contours.size() > 1)
        idxContour = findMatchingContour(contours,hierarchy,ptCentre);

    if (contours[idxContour].size() > 4)
    {
        cv::RotatedRect rectBody = fitEllipse(contours[idxContour]);
        // Make Slight Angle Corrections 0.5 Corrected
        iAngleOffset = rectBody.angle; //Make 0 Angle  the vertical image axis

        if  ((iAngleOffset) > 90) iAngleOffset-=180;
        if  ((iAngleOffset) < -90) iAngleOffset+=180;

        //Debug Show Shape Centre //
        //cv::circle(imgFishAnterior,rectBody.center,3,100,1);
        //Correct For Centre Offset

        ptCentreCorrection.x = ptCentre.x - rectBody.center.x;
        ptCentreCorrection.y = ptCentre.y - rectBody.center.y;
       // ptCentreCorrection.y //= tempCentre.x - rectBody.center.x;
        std::clog << "[info] scoreBlobReg: correct Template DAngle : " << iAngleOffset << " DX:" << ptCentreCorrection.x <<  std::endl;
    }


    ptRevised = ptCentreCorrection;
    if (abs(iAngleOffset-Angle) < 90 )
        RevisedAngle = (iAngleOffset+180)%360;
    else
        RevisedAngle = iAngleOffset;


    //cv::imshow("Fish Thresh",imgFishAnterior_thres);

    return iAngleOffset;

}

