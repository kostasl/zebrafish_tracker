///*
/// Uses Template Image To detect Matching and orientation of fish Body
/// \todo Need to Look into feature descriptors for matching as a faster method.
/// OpenCV includes some ready made, ORB being free. I Found an intuitive and external to OPENCV called FREAK
/// which could be implemented (uses a retina inspired sampling pattern to encode a binary string of the template)
///*



#include <template_detect.h>
#include <larvatrack.h>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
#include <opencv2/highgui/highgui.hpp>


/// CUDA //

#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    #include <opencv2/opencv_modules.hpp> //THe Cuda Defines are in here
    #include "opencv2/cudaimgproc.hpp"
    #include "opencv2/cudaarithm.hpp"
    #include <opencv2/core/cuda.hpp>
    #include <opencv2/photo/cuda.hpp>
    #include <opencv2/core/cuda_types.hpp>
#endif

#include <QString>
#include <random>
#include <QDirIterator>
#include <QDir>
#include <QDebug>





//#include <opencv2/cudaimgproc.hpp> //Template Matching \todo Compile OpenCv With CUDA Support

//extern double gTemplateMatchThreshold;
//extern int gFishTemplateAngleSteps;
//extern int gnumberOfTemplatesInCache;
extern cv::Mat gFishTemplateCache;
extern MainWindow* pwindow_main;
extern bool bTemplateSearchThroughRows;

#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    extern cv::Ptr<cv::cuda::TemplateMatching> gpu_MatchAlg;
#endif

static cv::Mat loadImage(const std::string& name)
{
    cv::Mat image = cv::imread(name,cv::IMREAD_GRAYSCALE );//
    if (image.empty())
    {
        std::cerr << "Can't load image - " << name << std::endl;
        exit(-1);
    }
    return image;
}


/// \brief Defines search area region and runs template matching
/// pt Cant be closer to image border than gFishBoundBoxSize - If it is it will fixed to this distance
/// Returns Angle of Matched Template, centre of Detected Template and Match Score
double doTemplateMatchAroundPoint(const cv::Mat& maskedImg_gray,cv::Point pt,int& iLastKnownGoodTemplateRow,int& iLastKnownGoodTemplateCol,
                                  int& detectedAngle,cv::Point& detectedPoint ,const cv::Mat& frameOut )
{
    /// Fix Bounds For Search point such that search temaplte region is not smaller than template size
    pt.x = (pt.x <= gTrackerState.gszTemplateImg.width)?(gTrackerState.gszTemplateImg.width/2): pt.x;
    pt.x = (maskedImg_gray.cols-pt.x <= gTrackerState.gszTemplateImg.width/2)?maskedImg_gray.cols-gTrackerState.gszTemplateImg.width/2: pt.x;

    pt.y = (pt.y <= gTrackerState.gszTemplateImg.height/2)?(gTrackerState.gszTemplateImg.height/2): pt.y;
    pt.y = (maskedImg_gray.rows-pt.y <= gTrackerState.gszTemplateImg.height/2)?maskedImg_gray.rows-gTrackerState.gszTemplateImg.height/2: pt.y;
    ///

    double maxMatchScore =0; //
    cv::Point gptmaxLoc; //point Of Bestr Match
    cv::Size szTempIcon(std::max(gTrackerState.gszTemplateImg.width,gTrackerState.gszTemplateImg.height),
                        std::max(gTrackerState.gszTemplateImg.width,gTrackerState.gszTemplateImg.height));
    assert(szTempIcon.width > 5 && szTempIcon.height> 5);
    //cv::Point rotCentre = cv::Point(gTrackerState.gszTemplateImg.height/2+4,gTrackerState.gszTemplateImg.width/2+2); //HACK for better positioning of anteriorFrame
    cv::Point rotCentre = cv::Point(szTempIcon.width/2,szTempIcon.height/2); //HACK for better positioning of anteriorFrame

    /// Check If Track Centre Point Contains An image that matches a fish template
    /// \todo make HeadPoint/Tail point a Propery of FishBlob
    //cv::Point centroid = fishblob->pt;
     //Locate Centroid Region at a point between blob Centroid And Detect HeadPoint on Curve
   // cv::Point centroid = ((cv::Point)fishblob->pt-gptHead)/3+gptHead;

    /// BOUND SEARCH REGION ///
    // Small Search Region When A Match has already been found
    cv::Point pBound1,pBound2;
    int iSearchRegionSize = max((int)(max(szTempIcon.height,szTempIcon.width)*0.1), (int)(0.2*gTrackerState.gFishBoundBoxSize));

    //Expand the Search Region If Fish Tracking Has been lost
    if (iLastKnownGoodTemplateRow == 0 && iLastKnownGoodTemplateCol == 0)
       iSearchRegionSize = 0.5*gTrackerState.gFishBoundBoxSize;


    pBound1 = cv::Point(std::max(0,std::min(maskedImg_gray.cols,pt.x-iSearchRegionSize)), std::max(0,std::min(maskedImg_gray.rows,pt.y-iSearchRegionSize)));
    pBound2 = cv::Point(std::max(0,std::min(maskedImg_gray.cols,pt.x+iSearchRegionSize)), std::max(0,std::min(maskedImg_gray.rows,pt.y+iSearchRegionSize)));

    // Look for Fish Template Within The Blob Region //
    cv::Rect rectFish(pBound1,pBound2);

    cv::Mat fishRegion(maskedImg_gray,rectFish); //Get Sub Region Image

    if (gTrackerState.bshowDetectorDebugImg)
        cv::imshow("template Fish region",fishRegion);

    //If blob exists but No Fish Model yet then Search Through Cache to improve matching;
    //bool findBestMatch = (vfishmodels.size() == 0);
    // Crude Way to detect if this fish has been searched before/ First Time this blob is searched?
    bool findBestMatch = (iLastKnownGoodTemplateRow == 0 && iLastKnownGoodTemplateCol == 0);
    if (findBestMatch)
        pwindow_main->LogEvent(QString("Search throughout templates for best match (slow)"));

   /// iLastKnownGoodTemplateRow will change to the row that matched the tracked larva
    int iLastTemplateRow = iLastKnownGoodTemplateRow;
    int AngleIdx = templatefindFishInImage(fishRegion,gFishTemplateCache,szTempIcon, maxMatchScore, gptmaxLoc,
                                           iLastKnownGoodTemplateRow,iLastKnownGoodTemplateCol,
                                           findBestMatch);

    //Log Change of Template Row
    if (iLastTemplateRow != iLastKnownGoodTemplateRow)
    {
        std::stringstream ss;
        ss << "Changed template row to: -> "  << iLastKnownGoodTemplateRow;
        pwindow_main->LogEvent(QString::fromStdString(ss.str()));
    }

    detectedAngle =AngleIdx*gTrackerState.gFishTemplateAngleSteps;

    //MaxLoc Coords Need inversion
    //gptmaxLoc.x = pBound2.x - gptmaxLoc.x;
    //gptmaxLoc.y = pBound2.y - gptmaxLoc.y;
    cv::Point top_left  = pBound1+gptmaxLoc; //Get top Left Corner Of Template Detected Region
    detectedPoint = top_left + rotCentre; //HACKED x1.5-Get Centre Of Template Detection Region - Used for Tracking
    /// Debug Draw //
    #ifdef _ZTFDEBUG_
        cv::rectangle(frameOut,rectFish,CV_RGB(200,0,0),1); //Ucomment to debug template search Region
        cv::circle(frameOut,top_left,3,CV_RGB(200,0,0),2); //Best Match Point in Region
    #endif



    return maxMatchScore;
}


///
/// \brief makeTemplateCache Creates a rotational replicate of the template across angles of size
/// stores them in consecutive blocks. Saves time from having to rotate the template EveryTime
///  Handles Symmetry Issues, If Template Is not Square, It copies to a square canvas before Rotation
/// \param templateIn
/// \param imgTemplateOut
/// \param iAngleStepDeg (minimum 1 degree)
///
void makeTemplateVar(cv::Mat& templateIn,cv::Mat& imgTemplateOut, int iAngleStepDeg)
{
    int iAngleIncrements = 360.0/iAngleStepDeg;
    int ifishtemplateAngle = 0;
    //Allocate Dark/Black Mat of appropriate Size to Fit All Replications of the templ.
    //Space Them Out Along Long Dim So no cutting Occurs
    int mxDim = std::max(templateIn.cols,templateIn.rows);

    assert(mxDim == max(gTrackerState.gszTemplateImg.height,gTrackerState.gszTemplateImg.width));

    imgTemplateOut = cv::Mat::zeros(mxDim,mxDim*iAngleIncrements,CV_8UC1);
    cv::Point tempCentre = cv::Point(templateIn.cols/2,templateIn.rows/2);
    cv::Point ptbottomRight = cv::Point(mxDim,mxDim); //Square fitting Largest Dimension
    cv::Point cntrdiff = cv::Point(mxDim/2,mxDim/2)-tempCentre ;
    cv::Rect templRegion(cv::Point(0,0),ptbottomRight); //Define A Rect region to snatch SubMatrices
    cv::Point rotCentre = cv::Point(templRegion.width/2,templRegion.height/2);

    //Find Alignment OffSet - Correct Template In Case Its Tilded Already //
    int iAngleOffset = 0;
    cv::Point ptCentreCorrection(0,0);
    /// Detect edges using Threshold
    cv::Mat templ_thres,template_blur;
    std::vector<std::vector<cv::Point> > contours;
    std::vector<cv::Vec4i> hierarchy;

    //Segment Fish Body Blob
    //cv::threshold( templateIn, templ_thres,thres, 255, cv::THRESH_BINARY );
    cv::GaussianBlur(templateIn,template_blur,cv::Size(7,7),5,5);
    cv::adaptiveThreshold(template_blur, templ_thres,255,cv::ADAPTIVE_THRESH_MEAN_C,cv::THRESH_BINARY,21,0); //Last Param Is const substracted from mean
    //ADAPTIVE_THRESH_MEAN_C

    /// Find contours
    cv::findContours( templ_thres, contours, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );

    //Should be one contour
    if (contours.size() == 1)
    {
        if (contours[0].size() > 4)
        {
            cv::RotatedRect rectBody = fitEllipse(contours[0]);
            // Make Slight Angle Corrections 0.5 Corrected
            iAngleOffset = rectBody.angle; //Make 0 Angle  the vertical image axis

            if  ((iAngleOffset) > 90) iAngleOffset-=180;
            if  ((iAngleOffset) < -90) iAngleOffset+=180;

            //Debug Show Shape Centre //
            cv::circle(templ_thres,rectBody.center,3,100,1);
            //Correct For Centre Offset

            ptCentreCorrection.x = tempCentre.x - rectBody.center.x;
           // ptCentreCorrection.y //= tempCentre.x - rectBody.center.x;
            std::clog << "[info] makeTemplateVar: correct Template DAngle : " << iAngleOffset << " DX:" << ptCentreCorrection.x <<  std::endl;
        }
    }else
    {
        std::clog << "[warning] makeTemplateVar: multiple contours detected on template " << std::endl;
    }


    // DEBUG //
    //std::string winname= QString(QString("TemplShape")+QString::number(gnumberOfTemplatesInCache)).toStdString();
    //cv::imshow(winname ,templ_thres);

    //Got through Each Angle
    for (int i=0;i<iAngleIncrements;i++)
    {
        //Make Rotation MAtrix
        cv::Mat Mrot = cv::getRotationMatrix2D(rotCentre,360.0-ifishtemplateAngle+iAngleOffset*0.5,1.0);

        //Copy to Larger Square

        //Get icon SubMatrix from Large Template Mat
        cv::Mat templ_rot = imgTemplateOut(templRegion);
        //ReLocate Template To Centre Of Large Canvas, Before Rotation
        cv::Point offset = cntrdiff+ptCentreCorrection;
        //Watch Out for Boundaries
        offset.x = (templateIn.cols+offset.x > templ_rot.cols)?0:max(0,offset.x);
        offset.y = (templateIn.rows+offset.y > templ_rot.rows)?0:max(0,offset.y);
        templateIn.copyTo(templ_rot(cv::Rect(0,0,templateIn.cols,templateIn.rows)+offset ));
        //Make Rotation Transformation
        cv::warpAffine(templ_rot,templ_rot,Mrot,templRegion.size());

        //Slide To Next Empty region box for the Next Angle
        templRegion.x +=templRegion.width;
        ifishtemplateAngle += iAngleStepDeg;
    }



}


///
/// \brief templatefindFishInImage Scans input image for templates and returns the best location
///  which that exceed a threshold value for a match
/// \param templRegion Rect of template img Size to look for within the larger template Cache
/// \param imgTemplCache Global large Mat 2D image having the array of Templates
/// \param templSz The size of each template icon as saved in the Cache - These are a square along largest templ Dim
/// \param startRow - Optimization So search begins from the most likely Template as of the last one
/// \param startCol - Optimization So search begins from the most likely Template Angle - Set To Zero And the All Angles Will be searched
/// \param findFirstMatch if true It Looks for 1st template row that exceeds threshold - otherwise it looks for best match through all cache
/// \note The calling Function needts reposition maxLoc To the global Frame, if imgGreyIn is a subspace of the image
/// if Row scanning is disabled when bTemplateSearchThroughRows is not set
/// Use of UMat for matchTemplate is superfluous , as the GPU is not Utilized and it causes a dealloc mem bug to trigger -  Removed UMat
/// A gpuAssisted  function for this is included in the bottom of the file.
int templatefindFishInImage(cv::Mat& imgRegionIn,cv::Mat& imgtemplCache,cv::Size templSz, double& matchScore,
                            cv::Point& ptBestMatchlocation,int& startRow,int& startCol,bool findFirstMatch)
{
  const int iIdxAngleMargin = 5; //Offset Of Angle To begin Before LastKnownGood Angle
  int matchColIdx = 0;
  int Colidx = 0; //Current Angle Index Being tested in the loop

  assert(!imgRegionIn.empty());
  assert(templSz.height*startRow <= imgtemplCache.rows);
  assert(templSz.width*startCol <= imgtemplCache.cols);

  //startRow = 0;
  cv::Mat templ_rot; //The Matched template
#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
        cv::cuda::GpuMat dimgRegionIn(imgRegionIn);
#endif


  int idRow = startRow;
  int ibestMatchRow  = startRow;
  double minVal;
  double maxVal = 0.0;
  double maxGVal = 0.0;
  cv::Point ptmaxLoc,ptminLoc;
  cv::Point ptGmaxLoc,ptGminLoc;
  cv::Mat outMatchConv; //Convolution Result
  //iAngles = imgtempl.cols/templRegion.width;
  //Slide over each little template icon


  cv::Point ptbottomRight = cv::Point(templSz.width,templSz.height);
  cv::Rect templRegion(cv::Point(0,0),ptbottomRight);

  // Start from Angle Region as set from last search
  if (startCol > iIdxAngleMargin)
      startCol -=iIdxAngleMargin; //Move to Template N Angle Steps anticlockwise

 //Best Match Flag Forces A Full Search Through the Cache
  if (findFirstMatch) //Reset Search To Start From Top
  {
      startRow = 0;
      idRow = 0;
      startCol = 0;
  }

  //Initialiaze At last known Good Location
  templRegion.x = std::max(0,std::min(templSz.width*startCol,imgtemplCache.cols));
  templRegion.y = std::max(0,std::min(templSz.height*startRow,imgtemplCache.rows));
  Colidx = startCol;

  ///Run Through All rotated Templates - optional starting row for optimization
  /// \note For Speed Up - We Remove Running Through All Rows - Stick to LastKnown Good Row- And Let Random switch when this fails take Care Of Changing Row
  int iScanRowLimit = imgtemplCache.rows;

  if (!gTrackerState.bTemplateSearchThroughRows) //Do not Search Subsequent Template Rows
      iScanRowLimit = std::min(templSz.height*startRow + templRegion.height,imgtemplCache.rows) ;

  // Run through croping/extracting each template sized window from the larger image
  for (int j=templSz.height*startRow; j<iScanRowLimit;j+=templRegion.height) //Remove for Speed Optimization.
  {
      templRegion.y    = j;
       /// Run Throught each  *Columns/Angle* (Ie Different Angles of this template
      //for (int i=templSz.width*startCol; i<imgtempl.cols;i+=templRegion.width)

      ///Limit Search To Within FIXED (15) Degrees/Templates from Starting Col.Point If Not Searching From The top
      while(templRegion.x < imgtemplCache.cols &&
            ((Colidx-startCol) < TEMPLATE_COL_SEARCH_REGION || startCol==0))
      {
        //Obtain next Template At Angle
         imgtemplCache(templRegion).copyTo(templ_rot) ;

        //Run Lib MatchTemplate COnvolution Either On CPU OR GPU //
         // WIth Current Setup Using GPU for multiple template Matching/scanning is too slow
#if defined(USE_CUDA) && defined(USE_CUDA_FOR_TEMPLATE_MATCHING) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
        maxVal = gpu_matchTemplate(templ_rot,dimgRegionIn,ptmaxLoc);
#else
         cv::matchTemplate(imgRegionIn,templ_rot,outMatchConv, cv::TM_CCORR_NORMED  ); // cv::TM_CCORR_NORMED CV_TM_CCOEFF_NORMED ,TM_SQDIFF_NORMED
         //Find Min Max Location

         //cv::flip(outMatchConv,outMatchConv,-1); //Flip H and V
         cv::minMaxLoc(outMatchConv,&minVal,&maxVal,&ptminLoc,&ptmaxLoc);
#endif         //Convolution  // CV_TM_SQDIFF_NORMED Poor Matching

        //Assume Value < 0.7 is non Fish,
       //maxVal   = 1-minVal; //When Using SQ Diff
       //ptmaxLoc = ptminLoc;
       if (maxGVal < maxVal)
       {
            maxGVal         = maxVal;
            ptGmaxLoc       =  cv::Point(ptmaxLoc.x, //+imgRegionIn.cols/2
                                         ptmaxLoc.y); //+imgRegionIn.rows/2 //The calling Function needs to reposition maxLoc To the global Frame
            matchColIdx     = Colidx;
            ibestMatchRow   = idRow;
            matchScore      = maxVal; //Save Score Of Best Match
            ptBestMatchlocation = ptGmaxLoc;

        }

        //Shift Region To Next Block
        if (findFirstMatch) //Skip Cols If Just Scanning
        {
            templRegion.x +=3*templSz.width;
            Colidx+=3;
        }
        else
        {
            templRegion.x +=templSz.width;
            Colidx++;
        }
      } //Loop Through Columns/Angle


        //Dont scan all Rows Just Check If Found on this One before Proceeding

       ///Check If Matching Exceeds threshold And Stop Loops - Return Found Column Idx//
       if (maxGVal >= gTrackerState.gTemplateMatchThreshold && maxGVal > 0.1) //!findFirstMatch
       {
           //Save Results To Output
           //matchScore    = maxGVal;
           //locations_tl  = ptGmaxLoc;
           startCol = matchColIdx; //Save AS Best - For next Iteration - Only If Exceeds Threshold
           break; ///Stop The loop Rows search Here

        }else //Reset StartCol Hint, So Next Time The whole Cache Row is scanned
            startCol = 0;

       idRow++;             //We Start Again From The Next Row
       Colidx               = 0;
       templRegion.x        = 0; //ReStart from 1st col


//   else{ //Nothing Found YEt-- Proceed To Next Template variation
//       matchColIdx  = 0;
//       matchScore   = maxGVal;
//       locations_tl = cv::Point(0,0);
//       //Didnt Find Template /Try Next Row
//       //startRow = 0;//Start From Top Of All Templates On Next Search
//   }


 } //Loop Through Rows - Starting From Last Call Best Match Row


  if (startRow != ibestMatchRow) //Store Best Row Match
  {
      std::stringstream ss;       // Log As Message //
      ss << "Ch. Templ. Row:" << startRow << " -> "  << ibestMatchRow;
      pwindow_main->LogEvent(QString::fromStdString(ss.str()));
      startRow = ibestMatchRow;
  }

// Check if template match passes user set threshold
 if (maxGVal < gTrackerState.gTemplateMatchThreshold)
 {
     gTrackerState.iTemplateMatchFailCounter++; //INcrease Count Of Failures
    // Choose next row Randomly
     startRow = (rand() % static_cast<int>(gTrackerState.gnumberOfTemplatesInCache - 0 + 1));//Start From RANDOM rOW On Next Search
     startCol = 0;

     std::stringstream ss;
//     // Log As Message //
     if (ibestMatchRow < gTrackerState.gnumberOfTemplatesInCache)
        ss << "Best row:" << ibestMatchRow << " but gives Low Match score:"<< maxGVal << ". Pick next Randomly ->"  << startRow;
     else
     {
         ss << "Reached End of Templates row:" << ibestMatchRow << " Start Over on next";
         startRow = 0;
     }
     pwindow_main->LogEvent(QString::fromStdString(ss.str()));

 }else{ // If template matched then stay on the same template row
     startRow = ibestMatchRow;

     /// ADJUST Threshold - Success - Move towards Match score
     if (gTrackerState.gTemplateMatchThreshold < gTrackerState.gTemplateMatchThreshold_UpLimit &&
             gTrackerState.gTemplateMatchThreshold >  gTrackerState.gTemplateMatchThreshold_LowLimit){
            //gTrackerState.gTemplateMatchThreshold += 0.01*(0.90*maxGVal-gTrackerState.gTemplateMatchThreshold);
            //gTrackerState.gTemplateMatchThreshold = std::min(gTrackerState.gTemplateMatchThreshold_UpLimit,
            //                                                 std::max(gTrackerState.gTemplateMatchThreshold_LowLimit, gTrackerState.gTemplateMatchThreshold));
            //pwindow_main->updateTemplateThres();
     }
     if (gTrackerState.bshowDetectorDebugImg)
        cv::imshow("TScore",outMatchConv);
     gTrackerState.iTemplateMatchFailCounter = 0; //Reset Counter of Failed Attempts

 }

 // Check if we are stuck with Too Many Template Match Fails, then Warn and lower match threshold.
 if (gTrackerState.iTemplateMatchFailCounter > gTrackerState.gnumberOfTemplatesInCache)
 {
     ///ADJUST THRESHOLD - Detect FAilures -
    if(gTrackerState.gTemplateMatchThreshold >  gTrackerState.gTemplateMatchThreshold_LowLimit)
    {
        //pwindow_main->LogEvent("[warning] Too many template match failures, lowering threshold.");
        //gTrackerState.gTemplateMatchThreshold +=  0.01*(0.90*maxGVal-gTrackerState.gTemplateMatchThreshold);
        //gTrackerState.gTemplateMatchThreshold = std::min(gTrackerState.gTemplateMatchThreshold_UpLimit,
        //                                                 std::max(gTrackerState.gTemplateMatchThreshold_LowLimit, gTrackerState.gTemplateMatchThreshold));

        pwindow_main->updateTemplateThres();
        //gTrackerState.iTemplateMatchFailCounter = 0; //Restart Counting With New Threshold
    }
 }

 if (templ_rot.cols > 0 && templ_rot.rows > 0)
         //cv::imshow("MTemplate",templ_rot);
        pwindow_main->showInsetTemplateimg(templ_rot); //Show The template Image Used
 else
     std::cerr << "template_detect: templ_rot is empty - matchcol:" << matchColIdx << " Row: " << ibestMatchRow << std::endl;

  return matchColIdx;
}


///
///\brief addTemplateToCache
///\note assumes all Templates are the same size
///
int addTemplateToCache(cv::Mat& imgTempl,cv::Mat& FishTemplateCache,int idxTempl)
{
    cv::Mat imgTempl_std(gTrackerState.gszTemplateImg.height,gTrackerState.gszTemplateImg.width,imgTempl.type());

    assert(imgTempl.rows <=  gTrackerState.gszTemplateImg.height &&
           imgTempl.cols <=  gTrackerState.gszTemplateImg.width);

    //Paste Template Into centre of Standard Template Frame Size
    cv::Rect pasteRegion((gTrackerState.gszTemplateImg.width-imgTempl.cols)/2,
                         (gTrackerState.gszTemplateImg.height-imgTempl.rows)/2
                         ,imgTempl.cols,imgTempl.rows);

    imgTempl.copyTo(imgTempl_std(pasteRegion));

    //Make Variations And store in template Cache
    cv::Mat fishTemplateVar;
    cv::Mat mtCacheRow,mtEnlargedCache;
    makeTemplateVar(imgTempl_std,fishTemplateVar, gTrackerState.gFishTemplateAngleSteps);

    ///Initialize The Cache if this the 1st Template added
    if (idxTempl == 0)
        FishTemplateCache = cv::Mat::zeros(fishTemplateVar.rows,fishTemplateVar.cols,CV_8UC1);
    else{
        //Copy COntents To Enlarged Cache and replace pointer
        mtEnlargedCache = cv::Mat::zeros(FishTemplateCache.rows+fishTemplateVar.rows,fishTemplateVar.cols,CV_8UC1);
        //Get Ref To Old Sized Cache Only
        mtCacheRow = mtEnlargedCache(cv::Rect(0,0,FishTemplateCache.cols,FishTemplateCache.rows));
        FishTemplateCache.copyTo(mtCacheRow); //Copy Old Cache into New replacing that part of empty cache
        mtEnlargedCache.copyTo(FishTemplateCache); //Copy Back So gFishTemplateCache = mtEnlargedCache;

        mtEnlargedCache.release();
        mtEnlargedCache.deallocate();
    }
     //Fill The Last (New Row) In The Cache
    mtCacheRow = FishTemplateCache(cv::Rect(0,fishTemplateVar.rows*(idxTempl),fishTemplateVar.cols,fishTemplateVar.rows));
    fishTemplateVar.copyTo(mtCacheRow); //Copy To Row In CAche
    gTrackerState.gnumberOfTemplatesInCache++; //Increment Count

    // DEBUG //
    //cv::imshow("Fish Template",FishTemplateCache(cv::Rect(0,0,std::max(imgTempl.cols,imgTempl.rows),FishTemplateCache.rows) ));
    //cv::imshow("Templ",imgTempl);
    //  //

    std::clog << "New Template added, Templ. Count now:" << gTrackerState.gnumberOfTemplatesInCache << std::endl;


    // Set Template Size
    //gTrackerState.gszTemplateImg.width = imgTempl.size().width; //Save TO Global Size Variable
    //gTrackerState.gszTemplateImg.height = imgTempl.size().height; //Save TO Global Size Variable

   return ++idxTempl;
}

///
/// \brief deleteLastTemplateRow
/// \param FishTemplateCache
/// \param idxTempl - Row To Remove
/// \return
///
int deleteTemplateRow(cv::Mat& imgTempl,cv::Mat& FishTemplateCache,int idxTempl)
{
    //Draw Black
    int mxDim = std::max(imgTempl.cols,imgTempl.rows);
    //\note RECT constructor takes starting point x,y, size_w, size_h)
    cv::Rect rectblankcv(0,mxDim*(idxTempl),FishTemplateCache.cols,mxDim);
    cv::Mat mFishTemplate_local;// = FishTemplateCache.getMat(cv::ACCESS_WRITE);
    cv::rectangle(mFishTemplate_local,rectblankcv,CV_RGB(0,0,0),cv::FILLED); //Blank It OUt

    if (idxTempl < 1)
    {
        pwindow_main->LogEvent("[Error] Attempted to delete Template with invalid ID (0)");
        return (0);
    }

    //If Removing Last Row, Then Its Simple
    //Shrink Template
    if (idxTempl == (gTrackerState.gnumberOfTemplatesInCache-1))
    {
        FishTemplateCache  = FishTemplateCache(cv::Rect(0,0,FishTemplateCache.cols,mxDim*(idxTempl)));
        gTrackerState.gnumberOfTemplatesInCache--;
    }else //Other Wise, We need to Cut and stich
    {
        //Cut In 2- Halves and rejoin
        cv::Mat mTop;
        cv::Mat mBottom;
        FishTemplateCache(cv::Rect(0,0,FishTemplateCache.cols,mxDim*(idxTempl))).copyTo(mTop) ;
        //cv::imshow("Template Top",mTop(cv::Rect(0,0,mxDim,mTop.rows)));

        FishTemplateCache(cv::Rect(0,mTop.rows+mxDim,mTop.cols,FishTemplateCache.rows-mTop.rows-mxDim)).copyTo(mBottom);
        //cv::imshow("Template Bottom",mBottom(cv::Rect(0,0,mxDim,mBottom.rows)));

        cv::Mat mtShrankCache   = cv::Mat::zeros(FishTemplateCache.rows-mxDim,FishTemplateCache.cols,CV_8UC1);

        gTrackerState.gnumberOfTemplatesInCache--;
        mTop.copyTo(mtShrankCache(cv::Rect(0,0,FishTemplateCache.cols,mxDim*(idxTempl))));
        mBottom.copyTo( mtShrankCache( cv::Rect(0,mTop.rows,FishTemplateCache.cols,mtShrankCache.rows-mTop.rows) ));
        mtShrankCache.copyTo(FishTemplateCache);
        mtShrankCache.deallocate();

    }

    // DEBUG //
    //cv::imshow("Updated Fish Template Cache",FishTemplateCache(cv::Rect(0,0,mxDim,FishTemplateCache.rows)));

    return 1;
}

/// Make Average Image For Classifier
cv::Mat makeMeanTemplateImage(std::vector<cv::Mat> vTemplImg)
{
    cv::Mat img_meanTempl = cv::Mat::zeros(vTemplImg[0].rows,vTemplImg[0].cols,CV_32FC(vTemplImg[0].channels()) ); //AccumWeight Result needs to be CV_32FC
    cv::Mat imgsample;
    std::vector<cv::Mat>::iterator it;
    for (it = vTemplImg.begin();it != vTemplImg.end();++it)
    {
        imgsample = (*it);
        try {
             cv::accumulateWeighted(imgsample*1.5,img_meanTempl,1.0/(float)vTemplImg.size());
        } catch (...) {
            pwindow_main->LogEvent("[Error] Making mean template image - check all templates match in size.");

        }

    }

    //img_meanTempl = img_meanTempl*(1.0/(float)vTemplImg.size());
    img_meanTempl.convertTo(img_meanTempl,vTemplImg[0].type());

    return (img_meanTempl);
}

/// \brief Loads all images from target dir and adds them to templateCache
std::vector<cv::Mat> loadTemplatesFromDirectory(QString strDir)
{
    QDir dirTempl(strDir);
    cv::Mat templFrame;
    std::vector<cv::Mat> vTempl_mat;

    int fileCount = 0;
    if (!dirTempl.exists())
    {
        qWarning("Cannot find the a template directory.");
        return vTempl_mat;
    }

    QStringList fileFilters; fileFilters << "*.png" << "*.tiff" << "*.pgm" << "*.png" << "*.jpg";
    QStringList imgFiles = QDir(strDir).entryList(fileFilters,QDir::Files,QDir::Name);
    strDir.append('/');
    QListIterator<QString> itfile (imgFiles);

    while (itfile.hasNext())
    {
      QString filename = itfile.next();
      std::string filepath = filename.prepend(strDir ).toStdString();

      qDebug() << "* Load Template: " << filename;
      templFrame  = loadImage(filepath);
      //Save to Glogal List
      vTempl_mat.push_back(templFrame);
      addTemplateToCache(templFrame,gFishTemplateCache,gTrackerState.gnumberOfTemplatesInCache);
      fileCount++;
    }


         qDebug() << "Loaded # " << fileCount << "Templates";

        return vTempl_mat;
}



//////////////////////// MATCH TEMPLATE EXAMPLE FOR GPU ////////////////////////
/// \brief gpu_matchTemplate
/// \param templ_h
/// \param image_h IMage region TO search In
/// \param ptBestMatch
/// \return MatchTemplateSCore
///
#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    double gpu_matchTemplate(cv::Mat templ_h,cv::cuda::GpuMat& dimage,cv::Point& ptBestMatch)
    {
    //cv::cuda::setDevice(0);	//initialize CUDA

    ////cv::Mat image_h = cv::imread(	"/home/buddy/Documents/workspace/OpenCVTemplateMatch1/src/input_image.jpg");
    ////cv::Mat templ_h = cv::imread(				"/home/buddy/Documents/workspace/OpenCVTemplateMatch1/src/template_image.jpg");

    cv::cuda::GpuMat dtempl(templ_h); //upload image on gpu

    cv::cuda::GpuMat dresult;

    //dtempl.upload(image_h);
    //cv::Ptr<cv::cuda::TemplateMatching> alg = cv::cuda::createTemplateMatching(templ_h.type(), CV_TM_CCORR_NORMED);


    //cv::cuda::GpuMat dst;
    gpu_MatchAlg->match(dimage, dtempl, dresult);

    //cv::cuda::normalize(dresult, dresult, 0, 1, cv::NORM_MINMAX, -1);
    double max_value;

    //cv::Point location;
    //Find Best Match
    cv::cuda::minMaxLoc(dresult, 0, &max_value, 0, &ptBestMatch);


    ////copying back to host memory for display
    //cv::Mat result_h;
    //result.download(result_h);

    //Return The Match Value
    return max_value;

    }
#endif
