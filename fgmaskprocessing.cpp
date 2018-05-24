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

#include <opencv2/core/ocl.hpp> //For setting setUseOpenCL

extern int g_Segthresh;
extern cv::Mat kernelOpen;
extern double dLearningRate; //Learning Rate During Initial BG Modelling done over MOGhistory frames
extern double dLearningRateNominal;
extern double gdMOGBGRatio;
//When Running Multiple Threads That Use BG Substractor - An SEGFault is hit in OpenCL
extern cv::Ptr<cv::BackgroundSubtractorMOG2> pMOG2; //MOG2 Background subtractor


extern ltROIlist vRoi;
extern cv::Point ptROI1;
extern cv::Point ptROI2; //This Default Value Is later Modified
extern cv::Size gszTemplateImg; //Used For ROI size

extern MainWindow* pwindow_main;

extern bool bStaticAccumulatedBGMaskRemove;

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
/// \brief updateBGFrame Update BG model for a fixed number of frames
/// \param frame
/// \param fgMask
/// \param nFrame
/// \return returns false when limit of updates is reached
///
bool updateBGFrame(cv::Mat& frame, cv::Mat& bgAcc, unsigned int nFrame,uint MOGhistory)
{


    std::vector<std::vector<cv::Point> > fishbodycontours;
    std::vector<cv::Vec4i> fishbodyhierarchy;
    bool ret = true;
    //Speed that stationary objects are removed
   // double dblRatioPxChanged    = 0.0;


    // Detect Food at Lower Thresh //
    cv::Mat bgMask,fgFishMask,fgFoodMask,frameImg_gray;

   // cv::equalizeHist( frame, frame );


    //Its Important to remove THe nOise Before doing MOG on the Pixels
    //cv::fastNlMeansDenoising(InputArray src, OutputArray dst, float h=3, int templateWindowSize=7, int searchWindowSize=21
///* Parameters:
/// src – Input 8-bit 1-channel, 2-channel or 3-channel image.
/// dst – Output image with the same size and type as src .
/// templateWindowSize – Size in pixels of the template patch that is used to compute weights. Should be odd. Recommended value 7 pixels
/// searchWindowSize – Size in pixels of the window that is used to compute weighted average for given pixel. Should be odd. Affect performance linearly: greater searchWindowsSize - greater denoising time. Recommended value 21 pixels
/// h – Parameter regulating filter strength. Big h value perfectly removes noise but also removes image details, smaller h value preserves details but also preserves some noise
////
    cv::fastNlMeansDenoising(frame, frameImg_gray,4.0,7, 41); /// \todo VS Vid 161 001 still fails in centre maybe increase the window size

    enhanceMask(frameImg_gray,bgMask,fgFishMask,fgFoodMask,fishbodycontours, fishbodyhierarchy);

    try
    {
        pMOG2->apply(frameImg_gray,fgFishMask,dLearningRate); //Let the Model Learn , Dont Interact With The Accumulated Mask
    }
    catch(...)
    {
        //##With OpenCL Support in OPENCV a Runtime Assertion Error Can occur /
        //In That case make OpenCV With No CUDA or OPENCL support
        //Ex: cmake -D CMAKE_BUILD_TYPE=RELEASE -D WITH_CUDA=OFF  -D WITH_OPENCL=OFF -D WITH_OPENCLAMDFFT=OFF -D WITH_OPENCLAMDBLAS=OFF -D CMAKE_INSTALL_PREFIX=/usr/local
        //A runtime Work Around Is given Here:
        std::clog << "MOG2 apply failed, probably multiple threads using OCL, switching OFF" << std::endl;
        pwindow_main->LogEvent("[Error] MOG2 failed, probably multiple threads using OCL, switching OFF");
        cv::ocl::setUseOpenCL(false); //When Running Multiple Threads That Use BG Substractor - An SEGFault is hit in OpenCL
    }


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

