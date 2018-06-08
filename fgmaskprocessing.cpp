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
#include <opencv2/core/cuda.hpp>
#include "opencv2/cudaimgproc.hpp"
#include "opencv2/cudaarithm.hpp"
#include <opencv2/photo/cuda.hpp>
#include <opencv2/core/cuda_types.hpp>


//const int gcFishContourSize         = ZTF_FISHCONTOURSIZE;

extern int g_Segthresh;
extern cv::Mat kernelOpen;
extern cv::Mat kernelDilateMOGMask;
extern cv::Mat kernelOpenfish;
extern cv::Point gptHead; ///\todo remove this global var hack

// Gaussian Curve Smoothing Kernels For fish Contour//
extern std::vector<double> gGaussian,dgGaussian,d2gGaussian; //These Are init. in main



extern const double dLearningRate; //Learning Rate During Initial BG Modelling done over MOGhistory frames
extern const double dLearningRateNominal;
extern double gdMOGBGRatio;
//When Running Multiple Threads That Use BG Substractor - An SEGFault is hit in OpenCL
extern cv::Ptr<cv::BackgroundSubtractorMOG2> pMOG2; //MOG2 Background subtractor

extern bool bRemovePixelNoise;
extern bool bUseBGModelling;
extern bool bUseGPU;
extern bool bUseOpenCL;
extern bool bshowMask;

extern ltROIlist vRoi;
extern cv::Point ptROI1;
extern cv::Point ptROI2; //This Default Value Is later Modified
extern cv::Size gszTemplateImg; //Used For ROI size

extern MainWindow* pwindow_main;

extern bool bStaticAccumulatedBGMaskRemove;

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

        cv::Mat frame;

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
        uint totFrames = capture.get(CV_CAP_PROP_FRAME_COUNT);
        //Do Not Overrun Vid Length In BG COmputation
        uint uiStopFrame = totFrames ; //std::min(totFrames,MOGhistoryLength);
        //read input data. ESC or 'q' for quitting
        uint uiLearnedFrames = 0;
        uint skipFrames = uiStopFrame/ MOGhistoryLength;

        window_main.setTotalFrames(uiStopFrame); //Set To BG Processing REgion

        while( !bExiting && (char)keyboard != 27 && nFrame < (uint) uiStopFrame && uiLearnedFrames < MOGhistoryLength)
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
            }
            else //Frame Grabbed - Process It
            {
                //Get Frame Position From Vid Sam
                nFrame = capture.get(CV_CAP_PROP_POS_FRAMES) + startFrameCount;
                window_main.nFrame = nFrame; //Update Window
                window_main.tickProgress();


                ///Make Global Roi on 1st frame if it doesn't prexist
                if (vRoi.size() == 0)
                {
                    ptROI2.x = frame.cols/2;
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


                if (bgAcc.empty()) //Make EMpty Mask
                    bgAcc = cv::Mat::zeros(frame.rows,frame.cols,CV_32FC(bgMask.channels()));

                frame.copyTo(frame,bgMask);
                cv::cvtColor( frame, frame, cv::COLOR_BGR2GRAY);

                updateBGFrame(frame, bgAcc, nFrame, MOGhistoryLength);
            }
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

           //Jump To Next Frame To Learn
           capture.set(CV_CAP_PROP_POS_FRAMES, nFrame+ skipFrames); //Move To Next So We Take MOGHistory Samples From the Video

        } //main While loop

        //Remove Low Values
        double uiMaxVal,uiMinVal;

        //Find Max Value,this should belong to stationary objects, and Use it as a relative measure to detect BG Objects
        cv::minMaxLoc(bgAcc,&uiMinVal,&uiMaxVal,0,0);
        cv::threshold(bgAcc,bgMask,uiMaxVal*0.05,255,cv::THRESH_BINARY); //All; Above 5% of Max are Stationary

        bgMask.convertTo(bgMask,CV_8UC1);

        //if (bStaticAccumulatedBGMaskRemove)
        //    cv::imshow("Accumulated BG Model",bgMask);
        pwindow_main->showVideoFrame(bgMask,nFrame);


        //delete capture object
        capture.release();




        //std::clog << gTimer.elapsed()/60000.0 << " Background Processing  loop. Finished" << std::endl;
         window_main.LogEvent(" Background Processing  loop. Finished");


      return nFrame;
} ///trackImageSequencefile




///
/// \brief processMasks Can Filter Pixel noise from frame_gray, Adds the BGModel Mask to a static mask, Can use Noise filtering to improve FG segmentation
/// \param frame_gray //Current greyScale Frame - Noise May be Removed If filtering Is Set To On
/// \param bgStaticMaskInOut The mask provided to processFrame, Includes Static Objects and ROI Region
///
void processMasks(cv::Mat& frame_gray,cv::Mat& bgStaticMaskInOut)
{
 cv::Mat fgMask;

#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
         if (bUseGPU)
         {
             dframe_gray.upload(frame_gray);

             if (bRemovePixelNoise)
             {        ///Remove Pixel Noise
                 ///* src – Input 8-bit 1-channel, 2-channel or 3-channel image.         ///        dst – Output image with the same size and type as src .         ///       templateWindowSize – Size in pixels of the template patch that is used to compute weights. Should be odd. Recommended value 7 pixels         ////        searchWindowSize – Size in pixels of the window that is used to compute weighted average for given pixel. Should be odd. Affect performance linearly: greater searchWindowsSize - greater denoising time. Recommended value 21 pixels         ///        h – Parameter regulating filter strength. Big h value perfectly removes noise but also removes image details, smaller h value preserves details but also preserves some noise         ///
                cv::cuda::fastNlMeansDenoising(dframe_gray, dframe_gray,2.0, 21,7);
                dframe_gray.download(frame_gray);
             }
             if (bUseBGModelling)
             {
                 try{
                        pMOG2->apply(dframe_gray,dframe_mask,dLearningRateNominal);
                    dframe_mask.download(fgMask);
                 }catch(...)
                 {
                     pwindow_main->LogEvent("[Error] CUDA MOG2 failed");
                     //cv::ocl::setUseOpenCL(false); //When Running Multiple Threads That Use BG Substractor - An SEGFault is hit in OpenCL
                 }
             } //BGMOdel
           }//Use GPU
         else{ //Note This Is the Same Code as In THe NO USE_CUDA Case
             if (bRemovePixelNoise)
                     cv::fastNlMeansDenoising(frame_gray, frame_gray,2.0,7, 21);
                   //Check If BG Ratio Changed
             if (bUseBGModelling)
             {
               try{
                   pMOG2->apply(frame_gray,fgMask,dLearningRateNominal);
               }catch(...)
               {
                   std::clog << "MOG2 apply failed, probably multiple threads using OCL, switching OFF" << std::endl;
                   pwindow_main->LogEvent("[Error] MOG2 failed, probably multiple threads using OCL, switching OFF");
                   cv::ocl::setUseOpenCL(false); //When Running Multiple Threads That Use BG Substractor - An SEGFault is hit in OpenCL
                   bUseOpenCL = false;
               }
             }//BGModel
         }
#else //NO GPU VERSION
      if (bRemovePixelNoise)
              cv::fastNlMeansDenoising(frame_gray, frame_gray,2.0,7, 21);
            //Check If BG Ratio Changed
      if (bUseBGModelling)
      {
        try{
            pMOG2->apply(frame_gray,fgMask,dLearningRateNominal);
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
      }//BGModel
#endif

      if (bUseBGModelling)
      {
        //Combine Masks and Remove Stationary Learned Pixels From Mask If Option Is Set
        if (bStaticAccumulatedBGMaskRemove && !bgStaticMaskInOut.empty() )//Although bgMask Init To zero, it may appear empty here!
        {
           //cv::bitwise_not(fgMask,fgMask);
            //gbMask Is Inverted Already So It Is The Accumulated FGMASK, and fgMask is the MOG Mask
           cv::bitwise_and(bgStaticMaskInOut,fgMask,fgMask); //Only On Non Stationary pixels - Ie Fish Dissapears At boundary
        }else
            fgMask.copyTo(bgStaticMaskInOut);
       }
      //NO FGMask As No Dynamic (MOG) Model Exists so simply Return the Static Mask


} //END PROCESSMASKS


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


///Remove Pixel Noise
///
///* src – Input 8-bit 1-channel, 2-channel or 3-channel image.
///        dst – Output image with the same size and type as src .
///       templateWindowSize – Size in pixels of the template patch that is used to compute weights. Should be odd. Recommended value 7 pixels
////        searchWindowSize – Size in pixels of the window that is used to compute weighted average for given pixel. Should be odd. Affect performance linearly: greater searchWindowsSize - greater denoising time. Recommended value 21 pixels
///        h – Parameter regulating filter strength. Big h value perfectly removes noise but also removes image details, smaller h value preserves details but also preserves some noise
///
//cv::fastNlMeansDenoising(InputArray src, OutputArray dst, float h=3, int templateWindowSize=7, int searchWindowSize=21

//frameImg_gray = frameImg.clone();
//cv::GaussianBlur(frameImg_gray,frameImg_blur,cv::Size(3,3),0);

outFishMask = cv::Mat::zeros(frameImg_gray.rows,frameImg_gray.cols,CV_8UC1);

//##CUDA For Dilate And Threshold
#if defined(USE_CUDA)
 if (bUseGPU)
 {
    cv::cuda::threshold(dframe_gray,dframe_thres,g_Segthresh,max_thresh,cv::THRESH_BINARY);
 }
 else
    cv::threshold( frameImg_gray, threshold_output, g_Segthresh, max_thresh, cv::THRESH_BINARY ); // Log Threshold Image + cv::THRESH_OTSU
#else
    cv::threshold( frameImg_gray, threshold_output, g_Segthresh, max_thresh, cv::THRESH_BINARY ); // Log Threshold Image + cv::THRESH_OTSU
#endif


//- Can Run Also Without THe BG Learning - But will detect imobile debri and noise MOG!
if (bUseBGModelling && !fgMask.empty()) //We Have a (MOG) Model In fgMask - So Remove those Stationary Pixels
{
    cv::Mat fgMask_dilate; //Expand The MOG Mask And Intersect with Threshold
    //cv::morphologyEx(fgMask,fgMask_dilate,cv::MORPH_OPEN,kernelOpenfish,cv::Point(-1,-1),1);
#if defined(USE_CUDA)
    if (bUseGPU) //Dilate The MOG Mask , and Combine
    {
        gpu_DilateFilter->apply(dframe_mask,dframe_mask); //Use Global GPU Mat to run dilation
        cv::cuda::bitwise_and(dframe_mask,dframe_thres,dframe_mask);
        dframe_mask.download(maskFGImg); //Transfer processed Mask Back To CPU Memory
    }else
    {
        cv::dilate(fgMask,fgMask_dilate,kernelDilateMOGMask,cv::Point(-1,-1),1);
        cv::bitwise_and(threshold_output,fgMask_dilate,maskFGImg); //Combine
    }
#else
    cv::dilate(fgMask,fgMask_dilate,kernelDilateMOGMask,cv::Point(-1,-1),1);
    cv::bitwise_and(threshold_output,fgMask_dilate,maskFGImg); //Combine
#endif
} //If BGModelling
else
{
    // Use thresh Only to Detect FG Fish Mask Detect
    //Returning The thresholded image is only required when No BGMask Exists
#if defined(USE_CUDA)
     if (bUseGPU)
        dframe_thres.download(threshold_output);
#endif

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

cv::morphologyEx(maskFGImg,outFoodMask,cv::MORPH_CLOSE,kernelOpen,cv::Point(-1,-1),1);

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
        if (bUseGPU)
            dframe_thres.download(threshold_output);

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
/// \brief updateBGFrame Update BG model for a fixed number of frames
/// \param frame
/// \param fgMask
/// \param nFrame
/// \return returns false when limit of updates is reached
///
bool updateBGFrame(cv::Mat& frameImg_gray, cv::Mat& bgAcc, unsigned int nFrame,uint MOGhistory)
{


    std::vector<std::vector<cv::Point> > fishbodycontours;
    std::vector<cv::Vec4i> fishbodyhierarchy;
    bool ret = true;
    //Speed that stationary objects are removed
   // double dblRatioPxChanged    = 0.0;


    // Detect Food at Lower Thresh //
    cv::Mat bgMask,fgFishMask,fgFoodMask;

   // cv::equalizeHist( frame, frame );


//#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)

//         dframe_gray.upload(frame);
//         cv::cuda::fastNlMeansDenoising(dframe_gray, dframe_gray,4.0, 41,7);
//         dframe_gray.download(frameImg_gray);
//#else
//    cv::fastNlMeansDenoising(frame, frameImg_gray,4.0,7, 41); /// \todo VS Vid 161 001 still fails in centre maybe increase the window size
//#endif

    processMasks(frameImg_gray,bgMask); //Applies MOG if bUseBGModelling is on
    enhanceMask(frameImg_gray,bgMask,fgFishMask,fgFoodMask,fishbodycontours, fishbodyhierarchy);

//    try
//    {
//        pMOG2->apply(frameImg_gray,fgFishMask,dLearningRate); //Let the Model Learn , Dont Interact With The Accumulated Mask
//    }
//    catch(...)
//    {
//        //##With OpenCL Support in OPENCV a Runtime Assertion Error Can occur /
//        //In That case make OpenCV With No CUDA or OPENCL support
//        //Ex: cmake -D CMAKE_BUILD_TYPE=RELEASE -D WITH_CUDA=OFF  -D WITH_OPENCL=OFF -D WITH_OPENCLAMDFFT=OFF -D WITH_OPENCLAMDBLAS=OFF -D CMAKE_INSTALL_PREFIX=/usr/local
//        //A runtime Work Around Is given Here:
//        std::clog << "MOG2 apply failed, probably multiple threads using OCL, switching OFF" << std::endl;
//        pwindow_main->LogEvent("[Error] MOG2 failed, probably multiple threads using OCL, switching OFF");
//        cv::ocl::setUseOpenCL(false); //When Running Multiple Threads That Use BG Substractor - An SEGFault is hit in OpenCL
//    }


    //Also Learn A pic of the stable features - Found In FoodMask - ie Fish Removed
    cv::accumulateWeighted(fgFoodMask,bgAcc,0.001);

    //dblRatioPxChanged = (double)cv::countNonZero(fgMask)/(double)fgMask.size().area();

    //DEBUG //
    //cv::imshow("fishMask",fgFishMask);

    pwindow_main->showVideoFrame(fgFishMask,nFrame);
    //cv::imshow("Accumulated Bg Model",bgAcc);

    //pMOG->apply(frame, fgMaskMOG,dLearningRate);
    //pGMG->apply(frame,fgMaskGMG,dLearningRate);


     //OPENCV 3 MORPHOLOGICAL
    //erode to get rid to food marks
    //cv::erode(fgMaskMOG2,fgMaskMOG2,kernel, cv::Point(-1,-1),3);
    //Do Close : erode(dilate())
    //cv::morphologyEx(fgMaskMOG2,fgMaskMOG2, cv::MORPH_CLOSE, kernelClose,cv::Point(-1,-1),2);
    //cv::dilate(fgMaskMOG2,fgMaskMOG2,kernel, cv::Point(-1,-1),4);
    //Apply Open Operation dilate(erode())
    //cv::morphologyEx(fgMaskMOG2,fgMaskMOG2, cv::MORPH_OPEN, kernel,cv::Point(-1,-1),2);



    return ret; //If False then tell calling function to stop updating
}

